/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "wfmngr.h"

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/processing/operator/transformation.h>
#include <seiscomp3/processing/operator/ncomps.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>
 
using namespace std;
using namespace Seiscomp;
using namespace Seiscomp::Processing;
using Seiscomp::Core::stringify;
using DataModel::ThreeComponents;

namespace {

template <class T>
T nextPowerOf2(T a, T min=1, T max=1<<31)
{
    int b = min;
    while (b < a)
    {
        b <<= 1;
        if (b > max)
            return -1;
    }
    return b;
}

}


namespace Seiscomp {
namespace HDD {

WfMngr::WfMngr(const std::string& recordStreamURL, const std::string& cacheDir, const std::string& wfDebugDir)
              : _recordStreamURL(recordStreamURL), _cacheDir(cacheDir) , _wfDebugDir(wfDebugDir)
{
    resetCounters();
}


std::string
WfMngr::getBandAndInstrumentCodes(const std::string& channelCode)
{
    if ( channelCode.size() >= 2)
        return channelCode.substr(0, 2);
    return "";
}


std::string
WfMngr::getOrientationCode(const std::string& channelCode)
{
    if ( channelCode.size() == 3)
        return channelCode.substr(2, 3);
    return "";
}


string 
WfMngr::waveformId(const Catalog::Phase& ph, const Core::TimeWindow& tw)
{
    return waveformId(ph.networkCode, ph.stationCode, ph.locationCode, ph.channelCode, tw);
}


string
WfMngr::waveformId(const string& networkCode, const string& stationCode,
                   const string& locationCode, const string& channelCode,
                   const Core::TimeWindow& tw)
{
    return stringify("%s.%s.%s.%s.%s.%s",
                     networkCode.c_str(), stationCode.c_str(),
                     locationCode.c_str(), channelCode.c_str(),
                     tw.startTime().iso().c_str(),
                     tw.endTime().iso().c_str());
}


Core::TimeWindow
WfMngr::traceTimeWindowToLoad(const Catalog::Phase& ph,
                              const Core::TimeWindow& neededTW,
                              bool useDiskCache,
                              bool performSnrCheck) const
{
    Core::TimeWindow twToLoad = neededTW;

    // if the SNR window is bigger than the xcorr window, than extend
    // the waveform time window 
    if ( performSnrCheck && _snr.minSnr > 0 )
    {
        Core::Time winStart = std::min( {
              neededTW.startTime(),
              ph.time + Core::TimeSpan(_snr.noiseStart),
              ph.time + Core::TimeSpan(_snr.signalStart)
        });
        Core::Time winEnd = std::max( {
                neededTW.endTime(),
                ph.time + Core::TimeSpan(_snr.noiseEnd),
                ph.time + Core::TimeSpan(_snr.signalEnd)
        });
        twToLoad = Core::TimeWindow(winStart, winEnd);
    }

    // Make sure to load at least 10 seconds of waveform. This avoid
    // re-loading waveforms for small changes in the settings, which
    // is frequent when the user is looking for the optimal configuration
    if ( useDiskCache )
    {
        const Core::TimeSpan additionalTime(5.);

        if ( twToLoad.startTime() > ph.time - additionalTime )
            twToLoad.setStartTime( ph.time - additionalTime );

        if ( twToLoad.endTime() < ph.time + additionalTime )
            twToLoad.setEndTime( ph.time + additionalTime );
    }

    return twToLoad;
}


/*
 * Return the waveform from the memory cache if present, otherwise load it
 */
GenericRecordCPtr
WfMngr::getWaveform(const Core::TimeWindow& tw,
                    const Catalog::Event& ev,
                    const Catalog::Phase& ph,
                    WfCache* memCache,
                    bool useDiskCache,
                    bool allowSnrCheck)
{
    string wfDesc = stringify("Waveform for Phase '%s' and Time slice from %s length %.2f sec",
                              string(ph).c_str(), tw.startTime().iso().c_str(), tw.length());

    bool doSnrCheck = (allowSnrCheck && _snr.minSnr > 0);

    const string wfId = WfMngr::waveformId(ph, tw);

    // try to load the waveform from the memory cache, if present
    if ( memCache )
    {
        // Check if the snr is good
        if ( ! doSnrCheck || _snrGoodWfs.count(wfId) != 0 )
        {
            const auto it = memCache->find(wfId);

            if ( it != memCache->end() )
            {
                return it->second; // waveform cached, just return it 
            }
        }
    }

    // Check if we have already excluded the trace because the snr is too high (save time)
    if ( doSnrCheck && _snrExcludedWfs.count(wfId) != 0 )
    {
        return nullptr;
    }

    // Check if we have already excluded the trace because we couldn't load it (save time)
    if ( _unloadableWfs.count(wfId) != 0 )
    {
        return nullptr;
    }

    //
    // Load the waveform, possibly perform a projection 123->ZNE or ZNE->ZRT,
    // filter it and finally save the result in the memory cache  for later re-use
    //
    bool projectionRequired = true;
    bool allComponents = false;
    DataModel::ThreeComponents tc;

    Core::Time refTime = tw.startTime();
    DataModel::SensorLocation *loc = Catalog::findSensorLocation(ph.networkCode, ph.stationCode, ph.locationCode, refTime);

    if ( ! loc )
    {
        // try to load waveform anyway, but no projection because we don't have the info
        projectionRequired = false;
    }
    else
    {
        string channelCodeRoot = WfMngr::getBandAndInstrumentCodes(ph.channelCode);
        allComponents = getThreeComponents(tc, loc, channelCodeRoot.c_str(), refTime);

        if ( ( tc.comps[ThreeComponents::Vertical] &&
               tc.comps[ThreeComponents::Vertical]->code() == ph.channelCode )        ||
             ( tc.comps[ThreeComponents::FirstHorizontal] &&
               tc.comps[ThreeComponents::FirstHorizontal]->code() == ph.channelCode ) ||
             ( tc.comps[ThreeComponents::SecondHorizontal] &&
               tc.comps[ThreeComponents::SecondHorizontal]->code() == ph.channelCode )
           )
        {
            projectionRequired = false;
        }
    }

    // if the SNR window is bigger than the xcorr window, than extend
    // the waveform time window
    const Core::TimeWindow twToLoad = traceTimeWindowToLoad(ph, tw, useDiskCache, doSnrCheck);

    // Load waveform:
    // - if no projection required, just load the requested component
    // - otherwise perform the projection 123->ZNE or ZNE->ZRT
    GenericRecordPtr trace;
    try {

        if ( ! projectionRequired )
        {
            trace = loadWaveform(twToLoad, ph.networkCode, ph.stationCode,
                                 ph.locationCode, ph.channelCode, useDiskCache);
        }
        else 
        {
            if ( ! allComponents )
            {
                SEISCOMP_DEBUG("Unable to fetch orientation information (%s)", wfDesc.c_str());
                _unloadableWfs.insert(wfId);
                _counters.wf_no_avail++;
                return nullptr;
            }
            trace = loadProjectWaveform(twToLoad, ev, ph, tc, loc, useDiskCache);
        }

    } catch ( exception &e ) {
        SEISCOMP_DEBUG("%s", e.what());
        _unloadableWfs.insert(wfId);
        _counters.wf_no_avail++;
        return nullptr;
    }

    // fitler waveform
    filter(*trace, true, _wfFilter.filterStr, _wfFilter.resampleFreq);

    // check SNR threshold
    if ( doSnrCheck  )
    {
        double snr = S2Nratio(trace, ph.time, _snr.noiseStart, _snr.noiseEnd, _snr.signalStart, _snr.signalEnd);
        if ( snr < _snr.minSnr ) 
        {
            _snrExcludedWfs.insert(wfId);
        }
        else
        {
            _snrGoodWfs.insert(wfId);
        }
    }

    // Trim waveform in case we loaded more data than requested (to compute SNR)
    if ( twToLoad != tw )
    {
        if ( ! WfMngr::trim(*trace, tw) )
        {
            SEISCOMP_DEBUG("Incomplete trace, not enough data (%s)", wfDesc.c_str());
            _unloadableWfs.insert(wfId);
            return nullptr;
        }
    }

    // save waveform into the cache
    if ( memCache )
    {
        (*memCache)[wfId] = trace;
    }

    // the trace has a high SNR, discard it if the SNR check was requested
    if ( doSnrCheck && _snrExcludedWfs.count(wfId) != 0 )
    {
        if ( _dump )
        {
            SEISCOMP_DEBUG("Trace has too low SNR, discard it (%s)", wfDesc.c_str());
            WfMngr::writeTrace(trace, waveformDebugPath(ev, ph, "snr-rejected") );
        }
        _counters.snr_low++;
        return nullptr;
    }

    if ( _dump )
    {
        string ext = ( ph.procInfo.source == Catalog::Phase::Source::THEORETICAL )
                   ? "theoretical" : (ph.isManual ? "manual" : "automatic" );
        WfMngr::writeTrace(trace, waveformDebugPath(ev, ph, ext) );
    }

    return trace; 
}


GenericRecordPtr
WfMngr::loadProjectWaveform(const Core::TimeWindow& tw,
                            const Catalog::Event& ev,
                            const Catalog::Phase& ph,
                            const DataModel::ThreeComponents& tc,
                            const DataModel::SensorLocation *loc,
                            bool useDiskCache) const
{
    string wfDesc = stringify("Waveform for Phase '%s' and Time slice from %s length %.2f sec",
                              string(ph).c_str(), tw.startTime().iso().c_str(), tw.length());

    SEISCOMP_DEBUG("Loading the 3 components waveforms (%s %s %s) to perform the projection...",
                   tc.comps[ThreeComponents::Vertical]->code().c_str(),
                   tc.comps[ThreeComponents::FirstHorizontal]->code().c_str(),
                   tc.comps[ThreeComponents::SecondHorizontal]->code().c_str());

    // orientation ZNE
    Math::Matrix3d orientationZNE;
    Math::Vector3d n;
    n.fromAngles(+deg2rad(tc.comps[ThreeComponents::Vertical]->azimuth()),
                 -deg2rad(tc.comps[ThreeComponents::Vertical]->dip())).normalize();
    orientationZNE.setColumn(2, n);
    n.fromAngles(+deg2rad(tc.comps[ThreeComponents::FirstHorizontal]->azimuth()),
                 -deg2rad(tc.comps[ThreeComponents::FirstHorizontal]->dip())).normalize();
    orientationZNE.setColumn(1, n);
    n.fromAngles(+deg2rad(tc.comps[ThreeComponents::SecondHorizontal]->azimuth()),
                 -deg2rad(tc.comps[ThreeComponents::SecondHorizontal]->dip())).normalize();
    orientationZNE.setColumn(0, n);

    // orientation ZRT
    Math::Matrix3d orientationZRT;
    double delta, az, baz;
    Math::Geo::delazi(ev.latitude, ev.longitude, loc->latitude(), loc->longitude(),
                      &delta, &az, &baz);
    orientationZRT.loadRotateZ(deg2rad(baz + 180.0));

    // transformation matrix
    Math::Matrix3d transformation;
    map<string,string> chCodeMap;

    string channelCodeRoot = getBandAndInstrumentCodes(ph.channelCode);
    string component = getOrientationCode(ph.channelCode);

    if ( component == "Z" || component == "N"  || component == "E" )
    {
        transformation = orientationZNE;
        chCodeMap[channelCodeRoot + "Z"] = tc.comps[ThreeComponents::Vertical]->code();
        chCodeMap[channelCodeRoot + "N"] = tc.comps[ThreeComponents::FirstHorizontal]->code();
        chCodeMap[channelCodeRoot + "E"] = tc.comps[ThreeComponents::SecondHorizontal]->code();

        SEISCOMP_DEBUG("Performing ZNE projection (channelCode %s -> %s) for %s",
            chCodeMap[ph.channelCode].c_str(), ph.channelCode.c_str(), wfDesc.c_str());
    }
    else if ( component == "R"  || component == "T" )
    {
        transformation.mult(orientationZRT, orientationZNE);
        //chCodeMap[channelCodeRoot + "Z"] = tc.comps[ThreeComponents::Vertical]->code();
        chCodeMap[channelCodeRoot + "R"] = tc.comps[ThreeComponents::FirstHorizontal]->code();
        chCodeMap[channelCodeRoot + "T"] = tc.comps[ThreeComponents::SecondHorizontal]->code();

        SEISCOMP_DEBUG("Performing ZRT projection (channelCode %s -> %s) for %s",
            chCodeMap[ph.channelCode].c_str(), ph.channelCode.c_str(), wfDesc.c_str()); 
    }
    else
    {
        string msg = stringify("Unknown channel '%s', cannot load %s", component.c_str(), wfDesc.c_str());
        throw runtime_error(msg); 
    }

    // Load the available components
    GenericRecordPtr tr1 = loadWaveform(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                                        tc.comps[ThreeComponents::Vertical]->code(), useDiskCache);
    GenericRecordPtr tr2 = loadWaveform(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                                        tc.comps[ThreeComponents::FirstHorizontal]->code(), useDiskCache);
    GenericRecordPtr tr3 = loadWaveform(tw, ph.networkCode, ph.stationCode, ph.locationCode,
                                        tc.comps[ThreeComponents::SecondHorizontal]->code(), useDiskCache);

    // The wrapper will direct 3 codes into the right slots using the
    // Stream configuration class and will finally use the transformation
    // operator. The advantage is that it will apply the configured gain for
    typedef Operator::StreamConfigWrapper<double,3,Operator::Transformation> OpWrapper;

    // Define the final operator class:
    //  1. Send channel codes to right slots
    //  2. Align 3 channels sample wise
    //  3. Transform the resulting 3 component trace with a rotation matrix
    typedef NCompsOperator<double,3,OpWrapper> Rotator;

    Processing::Stream streams[3];
    streams[2].init(tc.comps[ThreeComponents::Vertical]);
    streams[1].init(tc.comps[ThreeComponents::FirstHorizontal]);
    streams[0].init(tc.comps[ThreeComponents::SecondHorizontal]);
    Rotator op(OpWrapper(streams, Operator::Transformation<double,3>(transformation)));


    class DataStorer
    {
        public:
        DataStorer(const string& channelCode, const Core::TimeWindow& tw, map<string,string> chCodeMap)
        : chMap(chCodeMap[channelCode], channelCode)
        , _seq( new TimeWindowBuffer(tw) )
        {  }

        bool store(const Record *rec)
        {
            if (rec->channelCode() == chMap.first)
                _seq->feed(rec);
            return true;
        }

        const std::pair<string,string> chMap;
        std::shared_ptr<RecordSequence> _seq;
    };

    DataStorer projectedData(ph.channelCode, tw, chCodeMap);

    // The function that will be called after a transformed record was created
    //op.setStoreFunc(boost::bind(&RecordSequence::feed, seq, _1));
    op.setStoreFunc(boost::bind(&DataStorer::store, projectedData, _1));

    op.feed(tr1.get());
    op.feed(tr2.get());
    op.feed(tr3.get());

    std::shared_ptr<RecordSequence> seq = projectedData._seq;

    if ( seq->empty() )
    {
        string msg = stringify("No data after the projection for %s", wfDesc.c_str());
        throw runtime_error(msg);
    }

    GenericRecordPtr trace = new GenericRecord();

    if ( ! merge(*trace, *seq) )
    {
        string msg = stringify("Data records could not be merged into a single trace (%s)", wfDesc.c_str());
        throw runtime_error(msg);
    }

    trace->setChannelCode(ph.channelCode);

    if ( ! trim(*trace, tw) )
    {
        string msg = stringify("Incomplete trace, not enough data (%s)", wfDesc.c_str());
        throw runtime_error(msg);
    }

    return trace;
}


string
WfMngr::waveformDebugPath(const Catalog::Event& ev, const Catalog::Phase& ph, const std::string& ext) const
{
    string debugFile = stringify("ev%u.%s.%s.%s.%s.%s.%s.mseed", ev.id, 
                                 ph.networkCode.c_str(), ph.stationCode.c_str(),
                                 ph.locationCode.c_str(), ph.channelCode.c_str(),
                                 ph.type.c_str(), ext.c_str());
    return (boost::filesystem::path(_wfDebugDir)/debugFile).string();
}


string
WfMngr::waveformPath(const Catalog::Phase& ph, const Core::TimeWindow& tw) const
{
    return waveformPath(ph.networkCode, ph.stationCode, ph.locationCode, ph.channelCode, tw);
}


string
WfMngr::waveformPath(const string& networkCode, const string& stationCode,
                     const string& locationCode, const string& channelCode,
                     const Core::TimeWindow& tw) const
{
    string cacheFile = waveformId(networkCode, stationCode, locationCode, channelCode, tw) + ".mseed";
    return (boost::filesystem::path(_cacheDir)/cacheFile).string();
}


/*
 * Read a waveform from a chached copy on disk if present, otherwise
 * from the configured RecordStream
 */
GenericRecordPtr
WfMngr::loadWaveform(const Core::TimeWindow& tw,
                     const string& networkCode,
                     const string& stationCode,
                     const string& locationCode,
                     const string& channelCode,
                     bool useDiskCache) const
{
    string cacheFile = waveformPath(networkCode, stationCode, locationCode, channelCode, tw);

    GenericRecordPtr trace;
    // First try to read trace from disk cache
    if ( useDiskCache )
    {
        trace = readTrace(cacheFile);
        _counters.wf_cached++;
    }

    // if the trace is not cached then read it from the configured recordStream
    if ( !trace )
    {
        trace = readWaveformFromRecordStream(_recordStreamURL, tw, networkCode, stationCode, locationCode, channelCode);
        // then save the trace to disk for later usage
        if ( useDiskCache )
        {
            writeTrace(trace, cacheFile);
        }
        _counters.wf_downloaded++;
    }

    return trace;
}


GenericRecordPtr
WfMngr::readWaveformFromRecordStream(const string& recordStreamURL,
                                     const Core::TimeWindow& tw,
                                     const string& networkCode,
                                     const string& stationCode,
                                     const string& locationCode,
                                     const string& channelCode)
{
    IO::RecordStreamPtr rs = IO::RecordStream::Open( recordStreamURL.c_str() );
    if ( rs == nullptr )
    {
        string msg = "Cannot open RecordStream: " + recordStreamURL;
        throw runtime_error(msg);
    }

    rs->setTimeWindow(tw);
    rs->addStream(networkCode, stationCode, locationCode, channelCode);

    // Store each record in a RecordSequence
    IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
    std::shared_ptr<RecordSequence> seq( new TimeWindowBuffer(tw) );
    RecordPtr rec;
    while ( rec = inp.next() )
    {
        seq->feed(rec.get());
    }
    rs->close();

    if ( seq->empty() )
    {
        string msg = stringify("Data could not be loaded (stream %s.%s.%s.%s from %s length %.2f sec)",
                               networkCode.c_str(), stationCode.c_str(),
                               locationCode.c_str(), channelCode.c_str(),
                               tw.startTime().iso().c_str(), tw.length());
        throw runtime_error(msg);
    }

    GenericRecordPtr trace = new GenericRecord();

    if ( ! merge(*trace, *seq) )
    {
        string msg = stringify("Data records could not be merged into a single trace "
                               "(%s.%s.%s.%s from %s length %.2f sec)",
                               networkCode.c_str(), stationCode.c_str(),
                               locationCode.c_str(), channelCode.c_str(),
                               tw.startTime().iso().c_str(), tw.length());
        throw runtime_error(msg);
    }

    if ( ! trim(*trace, tw) )
    {
        string msg = stringify("Incomplete trace, not enough data for requested"
                               " time window (%s.%s.%s.%s from %s length %.2f sec)",
                               networkCode.c_str(), stationCode.c_str(),
                               locationCode.c_str(), channelCode.c_str(),
                               tw.startTime().iso().c_str(), tw.length());
        throw runtime_error(msg);
    }

    return trace;
}



bool
WfMngr::merge(GenericRecord &trace, const RecordSequence& seq)
{
    if ( seq.empty() )
    {
        return false;
    }

    RecordCPtr first = seq.front();
    RecordCPtr last;
    double samplingFrequency = first->samplingFrequency();
    Core::TimeSpan maxAllowedGap, maxAllowedOverlap;

    maxAllowedGap = Core::TimeSpan((double)(0.5 / samplingFrequency));
    maxAllowedOverlap = Core::TimeSpan((double)(-0.5 / samplingFrequency));

    trace.setNetworkCode(first->networkCode());
    trace.setStationCode(first->stationCode());
    trace.setLocationCode(first->locationCode());
    trace.setChannelCode(first->channelCode());

    trace.setStartTime(first->startTime());
    trace.setSamplingFrequency(samplingFrequency);

    Array::DataType datatype = first->data()->dataType();
    ArrayPtr arr = ArrayFactory::Create(datatype, datatype, 0, nullptr);

    for (const RecordCPtr& rec : seq )
    {
        if ( rec->samplingFrequency() != samplingFrequency ) {
            SEISCOMP_DEBUG("%s.%s.%s.%s: record sampling frequencies are not consistent: %f != %f",
                          trace.networkCode().c_str(),
                          trace.stationCode().c_str(),
                          trace.locationCode().c_str(),
                          trace.channelCode().c_str(),
                          samplingFrequency, rec->samplingFrequency());
            return false;
        }

        // Check for gaps and overlaps
        if ( last ) {
            Core::TimeSpan diff = rec->startTime()-last->endTime();
            if ( diff > maxAllowedGap ) {
                SEISCOMP_DEBUG("%s.%s.%s.%s: gap detected of %d.%06ds",
                              trace.networkCode().c_str(),
                              trace.stationCode().c_str(),
                              trace.locationCode().c_str(),
                              trace.channelCode().c_str(),
                              (int)diff.seconds(), (int)diff.microseconds());
                return false;
            }

            if ( diff < maxAllowedOverlap ) {
                SEISCOMP_DEBUG("%s.%s.%s.%s: overlap detected of %fs",
                              trace.networkCode().c_str(),
                              trace.stationCode().c_str(),
                              trace.locationCode().c_str(),
                              trace.channelCode().c_str(),
                              (double)diff);
                return false;
            }
        }

        arr->append( (Array *)(rec->data()));

        last = rec;
    }

    trace.setData(arr.get());

    return true;
}


bool WfMngr::trim(GenericRecord &trace, const Core::TimeWindow& tw)
{
    int ofs = (int)(double(tw.startTime() - trace.startTime())*trace.samplingFrequency());
    int samples = (int)(tw.length()*trace.samplingFrequency());

    // Not enough data at start of time window
    if ( ofs < 0 )
    {
        SEISCOMP_DEBUG("%s: need %d more samples in past",
                       trace.streamID().c_str(), -ofs);
        return false;
    }

    // Not enough data at end of time window
    if ( ofs+samples > trace.data()->size() )
    {
        SEISCOMP_DEBUG("%s: need %d more samples past the end",
                       trace.streamID().c_str(), trace.data()->size()-samples-ofs);
        return false;
    }

    ArrayPtr sliced = trace.data()->slice(ofs, ofs+samples);

    trace.setStartTime(tw.startTime());
    trace.setData(sliced.get());

    return true;
}



void WfMngr::filter(GenericRecord &trace, bool demeaning, const std::string& filterStr, double resampleFreq)
{
    DoubleArray *data = DoubleArray::Cast(trace.data());

    if (demeaning)
    {
        *data -= data->mean();
        trace.dataUpdated();
    }

    if ( resampleFreq > 0)
    {
        resample(trace, resampleFreq, true);
    }

    if ( ! filterStr.empty() )
    {
        string filterError;
        auto filter = Math::Filtering::InPlaceFilter<double>::Create(filterStr, &filterError);
        if ( !filter )
        {
            string msg = stringify("Filter creation failed %s: %s", filterStr.c_str(), filterError.c_str());
            throw runtime_error(msg);
        }
        filter->setSamplingFrequency(trace.samplingFrequency());
        filter->apply(data->size(), data->typedData());
        delete filter;
        trace.dataUpdated();
    }
}


void WfMngr::resample(GenericRecord &trace, double sf, bool average)
{
    if ( sf <= 0 )
        return;

    if ( trace.samplingFrequency() == sf )
        return;

    DoubleArray *data = DoubleArray::Cast(trace.data());
    double step = trace.samplingFrequency() / sf;

    if ( trace.samplingFrequency() < sf ) // upsampling
    {
        double fi = data->size() - 1;
        data->resize( data->size() / step );

        for( int i = data->size()-1; i >= 0; i-- )
        {
            (*data)[i] = (*data)[(int)fi];
            fi -= step;
        }
    }
    else // downsampling
    {
        int w = average?step*0.5 + 0.5:0;
        int i = 0;
        double fi = 0.0;
        int cnt = data->size();

        if ( w <= 0 )
        {
            while ( fi < cnt ) {
                (*data)[i++] = (*data)[(int)fi];
                fi += step;
            }
        }
        else
        {
            while ( fi < cnt )
            {
                int ci = (int)fi;
                double scale = 1.0;
                double v = (*data)[ci];

                for ( int g = 1; g < w; ++g )
                {
                    if ( ci >= g )
                    {
                        v += (*data)[ci-g];
                        scale += 1.0;
                    }

                    if ( ci+g < cnt )
                    {
                        v += (*data)[ci+g];
                        scale += 1.0;
                    }
                }

                v /= scale;

                (*data)[i++] = v;
                fi += step;
            }
        }
        data->resize(i);
    }
    trace.setSamplingFrequency((double)sf);
    trace.dataUpdated();
}


double WfMngr::S2Nratio(const GenericRecordCPtr& tr, const Core::Time& pickTime,
                 double noiseOffsetStart, double noiseOffsetEnd,
                 double signalOffsetStart, double signalOffsetEnd)
{
    const double *data = DoubleArray::ConstCast(tr->data())->typedData();
    const int size = tr->data()->size();
    const double freq = tr->samplingFrequency();
    const Core::Time dataStartTime = tr->startTime();

    // convert time w.r.t. guiding pick time to sample number
    auto secToSample = [&, freq, size](double sec) { 
        return std::min(std::max(std::round(sec * freq), 0.), size-1.);
    };
    const double pickOffset = (pickTime - dataStartTime).length();
    const int noiseStart  = secToSample(noiseOffsetStart  + pickOffset);
    const int noiseEnd    = secToSample(noiseOffsetEnd    + pickOffset);
    const int signalStart = secToSample(signalOffsetStart + pickOffset);
    const int signalEnd   = secToSample(signalOffsetEnd   + pickOffset);

    if ( (std::min({noiseStart,noiseEnd,signalStart,signalEnd}) < 0)    ||
         (std::max({noiseStart,noiseEnd,signalStart,signalEnd}) >= size) )
    {
        SEISCOMP_ERROR("Cannot compute S2N ratio: noise/signal windows exceed waveform boundaries");
        return -1;
    }

    // Get maximum (absolute) amplitude in noise window:
    double noiseMax = -1.0;
    for (int i = noiseStart; i < noiseEnd; i++)
    {
        noiseMax = std::max(std::abs(data[i]), noiseMax);
    }

    // Get maximum (absolute) amplitude in signal window:
    double signalMax = -1.0;
    for (int i = signalStart; i < signalEnd; i++)
    {
        signalMax = std::max(std::abs(data[i]), signalMax);
    }

    return signalMax/noiseMax;
}


void WfMngr::writeTrace(GenericRecordCPtr trace, const std::string& file)
{
    if ( ! trace )
        return;

    try {
        std::ofstream ofs(file);
        IO::MSeedRecord msRec(*trace);
        int reclen = msRec.data()->size()*msRec.data()->bytes() + 64;
        reclen = nextPowerOf2<int>(reclen, 128, 1048576); // MINRECLEN 128, MAXRECLEN 1048576
        if (reclen > 0)
        {
            msRec.setOutputRecordLength(reclen);
            msRec.write(ofs);
        }
    } catch ( exception &e ) {
        SEISCOMP_WARNING("Couldn't write waveform to disk %s: %s", file.c_str(), e.what());
    }
}


GenericRecordPtr WfMngr::readTrace(const std::string& file)
{
    if ( ! Util::fileExists(file) )
        return nullptr;

    try {
        std::ifstream ifs(file);
        IO::MSeedRecord msRec(Array::DOUBLE, Record::Hint::DATA_ONLY);
        msRec.read(ifs);
        GenericRecordPtr trace = new GenericRecord(msRec);
        trace->setData(msRec.data()->clone()); // copy data too
        return trace;
    } catch ( exception &e ) {
        SEISCOMP_WARNING("Couldn't load waveform %s: %s", file.c_str(), e.what());
        return nullptr;
    }
}


} // HDD
} // Seiscomp
