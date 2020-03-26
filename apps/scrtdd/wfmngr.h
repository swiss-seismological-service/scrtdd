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

#ifndef __RTDD_APPLICATIONS_WFMNGR_H__
#define __RTDD_APPLICATIONS_WFMNGR_H__

#include "catalog.h"

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/datamodel/utils.h>

#include <set>
#include <map>
#include <vector>

namespace Seiscomp {
namespace HDD {


DEFINE_SMARTPOINTER(WfMngr);

class WfMngr : public Core::BaseObject {

    public:

        WfMngr(const std::string& recordStreamURL, const std::string& cacheDir, const std::string& wfDebugDir);
        virtual ~WfMngr() { }

        void setWaveformDebug(bool dump ) { _dump = dump; }

        void setSnr(double minSnr, double noiseStart, double noiseEnd, double signalStart, double signalEnd)
        {
            _snr.minSnr      = minSnr;
            _snr.noiseStart  = noiseStart;
            _snr.noiseEnd    = noiseEnd;
            _snr.signalStart = signalStart;
            _snr.signalEnd   = signalEnd;
        }

        void setProcessing(const std::string& filterStr, double resampleFreq)
        {
            _wfFilter.filterStr    = filterStr;
            _wfFilter.resampleFreq = resampleFreq;
        }

        void resetCounters() { _counters = {0}; }

        void getCounters(unsigned& snr_low, unsigned& wf_no_avail, unsigned& wf_cached, unsigned& wf_downloaded)
        { 
            snr_low       = _counters.snr_low;
            wf_no_avail   = _counters.wf_no_avail;
            wf_cached     = _counters.wf_cached;
            wf_downloaded = _counters.wf_downloaded;
        }

        typedef std::map<std::string, GenericRecordCPtr> WfCache;

        GenericRecordCPtr getWaveform(const Core::TimeWindow& tw,
                                     const Catalog::Event& ev,
                                     const Catalog::Phase& ph,
                                     std::map<std::string,GenericRecordCPtr>* memCache,
                                     bool useDiskCache,
                                     bool allowSnrCheck);
        //
        //  static: utility functions
        //
        static GenericRecordPtr readWaveformFromRecordStream(const std::string& recordStreamURL,
                                                      const Core::TimeWindow& tw,
                                                      const std::string& networkCode,
                                                      const std::string& stationCode,
                                                      const std::string& locationCode,
                                                      const std::string& channelCode);

        static bool merge(GenericRecord &trace, const RecordSequence& seq);
        static bool trim(GenericRecord &trace, const Core::TimeWindow& tw);
        static void filter(GenericRecord &trace, bool demeaning=true, const std::string& filterStr="", double resampleFreq=0);
        static void resample(GenericRecord& trace, double sf, bool average);

        static void writeTrace(GenericRecordCPtr trace, const std::string& file);
        static GenericRecordPtr readTrace(const std::string& file);

        static std::string getBandAndInstrumentCodes(const std::string& channelCode);
        static std::string getOrientationCode(const std::string& channelCode);

        static double S2Nratio(const GenericRecordCPtr& tr, const Core::Time& guidingPickTime,
                              double noiseOffsetStart, double noiseOffsetEnd,
                              double signalOffsetStart, double signalOffsetEnd);
    private:

        GenericRecordPtr loadWaveform(const Core::TimeWindow& tw,
                                      const std::string& networkCode,
                                      const std::string& stationCode,
                                      const std::string& locationCode,
                                      const std::string& channelCode,
                                      bool useDiskCache) const;
        GenericRecordPtr loadProjectWaveform(const Core::TimeWindow& tw,
                                             const Catalog::Event& ev,
                                             const Catalog::Phase& ph,
                                             const DataModel::ThreeComponents& tc,
                                             const DataModel::SensorLocation *loc,
                                             bool useDiskCache) const;

        std::string waveformPath(const Catalog::Phase& ph, const Core::TimeWindow& tw) const;

        std::string waveformPath(const std::string& networkCode, const std::string& stationCode,
                                 const std::string& locationCode, const std::string& channelCode,
                                 const Core::TimeWindow& tw) const;

        std::string waveformDebugPath(const Catalog::Event& ev, const Catalog::Phase& ph,
                                      const std::string& ext) const;

        Core::TimeWindow traceTimeWindowToLoad(const Catalog::Phase& ph,
                                               const Core::TimeWindow& neededTW,
                                               bool useDiskCache,
                                               bool performSnrCheck) const;

        static std::string waveformId(const Catalog::Phase& ph, const Core::TimeWindow& tw);
        static std::string waveformId(const std::string& networkCode, const std::string& stationCode,
                                     const std::string& locationCode, const std::string& channelCode,
                                     const Core::TimeWindow& tw);

    private:
        std::string _recordStreamURL;
        std::string _cacheDir;
        std::string _wfDebugDir;
        std::set<std::string> _unloadableWfs;
        std::set<std::string> _snrGoodWfs;
        std::set<std::string> _snrExcludedWfs;

        bool _dump = false;

        struct {
            double minSnr = 0;
            double noiseStart = 0;
            double noiseEnd = 0;
            double signalStart = 0;
            double signalEnd = 0;
        } _snr;

        struct {
            std::string filterStr = "";
            double resampleFreq = 0;
        } _wfFilter;

        struct {
            unsigned snr_low;
            unsigned wf_no_avail;
            unsigned wf_cached;
            unsigned wf_downloaded;
        } mutable _counters;
};

}
}

#endif
