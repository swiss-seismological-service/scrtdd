/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as          *
 * published by the Free Software Foundation, either version 3 of the      *
 * License, or (at your option) any later version.                         *
 *                                                                         *
 * This software is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __HDD_WAVEFORM_H__
#define __HDD_WAVEFORM_H__

#include "catalog.h"

#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/datamodel/utils.h>
#include <seiscomp3/io/recordstream.h>

#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace Seiscomp {
namespace HDD {
namespace Waveform {

GenericRecordPtr
readWaveformFromRecordStream(const std::string &recordStreamURL,
                             const Core::TimeWindow &tw,
                             const std::string &networkCode,
                             const std::string &stationCode,
                             const std::string &locationCode,
                             const std::string &channelCode,
                             double tolerance,
                             double minAvailability);

bool projectionRequired(const Core::TimeWindow &tw,
                        const Catalog::Phase &ph,
                        const Catalog::Event &ev,
                        DataModel::ThreeComponents &tc,
                        DataModel::SensorLocation *&loc);

GenericRecordPtr projectWaveform(const Core::TimeWindow &tw,
                                 const Catalog::Phase &ph,
                                 const Catalog::Event &ev,
                                 double tolerance,
                                 double minAvailability,
                                 const GenericRecordCPtr &tr1,
                                 const GenericRecordCPtr &tr2,
                                 const GenericRecordCPtr &tr3,
                                 const DataModel::ThreeComponents &tc,
                                 const DataModel::SensorLocation *loc);

GenericRecordPtr contiguousRecord(const RecordSequence &seq,
                                  const Core::TimeWindow &tw,
                                  double tolerance,
                                  double minAvailability);

bool trim(GenericRecord &trace, const Core::TimeWindow &tw);

void filter(GenericRecord &trace,
            bool demeaning               = true,
            const std::string &filterStr = "",
            double resampleFreq          = 0);
void resample(GenericRecord &trace, double new_sf);

void writeTrace(GenericRecordCPtr trace, const std::string &file);
GenericRecordPtr readTrace(const std::string &file);

double computeSnr(const GenericRecordCPtr &tr,
                  const Core::Time &pickTime,
                  double noiseOffsetStart,
                  double noiseOffsetEnd,
                  double signalOffsetStart,
                  double signalOffsetEnd);

bool xcorr(const GenericRecordCPtr &tr1,
           const GenericRecordCPtr &tr2,
           double maxDelay,
           bool qualityCheck,
           double &delayOut,
           double &coeffOut);

void crossCorrelation(const double *dataS,
                      const int sizeS,
                      const double *dataL,
                      const int sizeL,
                      bool qualityCheck,
                      double &delayOut,
                      double &coeffOut);

std::string getBandAndInstrumentCodes(const std::string &channelCode);
std::string getOrientationCode(const std::string &channelCode);

std::string waveformId(const Core::TimeWindow &tw,
                       const std::string &networkCode,
                       const std::string &stationCode,
                       const std::string &locationCode,
                       const std::string &channelCode);
std::string waveformId(const HDD::Catalog::Phase &ph,
                       const Core::TimeWindow &tw);

DEFINE_SMARTPOINTER(Loader);

class Loader : public Core::BaseObject
{

public:
  Loader(const std::string &recordStream) : _recordStreamURL(recordStream) {}

  virtual ~Loader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev);

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev,
                                bool demeaning,
                                const std::string &filterStr,
                                double resampleFreq);

  // empty string disables debugging
  void setDebugDirectory(const std::string &directory = "")
  {
    _wfDebugDir = directory;
  }

  unsigned _counters_wf_no_avail   = 0;
  unsigned _counters_wf_downloaded = 0;

protected:
  GenericRecordPtr process(const GenericRecordCPtr &trace,
                           bool demeaning,
                           const std::string &filterStr,
                           double resampleFreq);

  const std::string _recordStreamURL;
  std::string _wfDebugDir;

  static constexpr double _tolerance       = 0.1;
  static constexpr double _minAvailability = 0.95;
};

class CompositeLoader : public Loader
{

public:
  CompositeLoader(LoaderPtr auxLdr) : Loader("") { setAuxLoader(auxLdr); }

  virtual ~CompositeLoader() {}

  void setAuxLoader(LoaderPtr auxLdr)
  {
    if (!auxLdr) throw std::runtime_error("Auxiliary loader cannot be null");
    _auxLdr = auxLdr;
  }

protected:
  LoaderPtr _auxLdr;
};

DEFINE_SMARTPOINTER(DiskCachedLoader);

class DiskCachedLoader : public CompositeLoader
{
public:
  DiskCachedLoader(LoaderPtr auxLdr, const std::string &cacheDir)
      : CompositeLoader(auxLdr), _cacheDir(cacheDir)
  {}

  virtual ~DiskCachedLoader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev);

  bool isCached(const Core::TimeWindow &tw,
                const Catalog::Phase &ph,
                const Catalog::Event &ev);

  unsigned _counters_wf_cached = 0;

protected:
  virtual GenericRecordCPtr getFromCache(const Core::TimeWindow &tw,
                                         const std::string &networkCode,
                                         const std::string &stationCode,
                                         const std::string &locationCode,
                                         const std::string &channelCode);

  virtual void storeInCache(const Core::TimeWindow &tw,
                            const std::string &networkCode,
                            const std::string &stationCode,
                            const std::string &locationCode,
                            const std::string &channelCode,
                            const GenericRecordCPtr &trace);

  std::string waveformPath(const std::string &cacheDir,
                           const Core::TimeWindow &tw,
                           const std::string &networkCode,
                           const std::string &stationCode,
                           const std::string &locationCode,
                           const std::string &channelCode) const;

  std::string _cacheDir;
};

DEFINE_SMARTPOINTER(MemCachedLoader);

class MemCachedLoader : public CompositeLoader
{
public:
  MemCachedLoader(LoaderPtr auxLdr) : CompositeLoader(auxLdr) {}

  virtual ~MemCachedLoader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev)
  {
    throw std::runtime_error("Cannot return unprocessed data");
  }

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev,
                                bool demeaning,
                                const std::string &filterStr,
                                double resampleFreq);

  bool isCached(const Core::TimeWindow &tw,
                const Catalog::Phase &ph,
                const Catalog::Event &ev);

  unsigned _counters_wf_cached = 0;

protected:
  virtual GenericRecordCPtr getFromCache(const Core::TimeWindow &tw,
                                         const std::string &networkCode,
                                         const std::string &stationCode,
                                         const std::string &locationCode,
                                         const std::string &channelCode);

  virtual void storeInCache(const Core::TimeWindow &tw,
                            const std::string &networkCode,
                            const std::string &stationCode,
                            const std::string &locationCode,
                            const std::string &channelCode,
                            const GenericRecordCPtr &trace);

  std::unordered_map<std::string, GenericRecordCPtr> _waveforms;
};

DEFINE_SMARTPOINTER(ExtraLenLoader);

class ExtraLenLoader : public CompositeLoader
{
public:
  ExtraLenLoader(LoaderPtr auxLdr, double traceMinLen)
      : ExtraLenLoader(auxLdr, traceMinLen / 2, traceMinLen / 2)
  {}

  ExtraLenLoader(LoaderPtr auxLdr, double beforePickLen, double afterPickLen)
      : CompositeLoader(auxLdr), _beforePickLen(beforePickLen),
        _afterPickLen(afterPickLen)
  {}

  virtual ~ExtraLenLoader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev);

  Core::TimeWindow traceTimeWindowToLoad(const Core::TimeWindow &neededTW,
                                         const Core::Time &pickTime) const;

protected:
  double _beforePickLen; // secs
  double _afterPickLen;  // secs
};

DEFINE_SMARTPOINTER(SnrFilteredLoader);

class SnrFilteredLoader : public CompositeLoader
{

public:
  SnrFilteredLoader(LoaderPtr auxLdr,
                    double minSnr,
                    double noiseStart,
                    double noiseEnd,
                    double signalStart,
                    double signalEnd)
      : CompositeLoader(auxLdr), _snr{minSnr, noiseStart, noiseEnd, signalStart,
                                      signalEnd}
  {}

  virtual ~SnrFilteredLoader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev)
  {
    throw std::runtime_error("Cannot compute SNR on unprocessed data");
  }

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev,
                                bool demeaning,
                                const std::string &filterStr,
                                double resampleFreq);

  Core::TimeWindow snrTimeWindow(const Core::Time &pickTime) const;

  bool goodSnr(const GenericRecordCPtr &trace,
               const Core::Time &pickTime) const;

  unsigned _counters_wf_snr_low = 0;

protected:
  struct
  {
    double minSnr;
    double noiseStart;  // secs relative to pick time
    double noiseEnd;    // secs relative to pick time
    double signalStart; // secs relative to pick time
    double signalEnd;   // secs relative to pick time
  } _snr;

  std::unordered_set<std::string> _snrGoodWfs;
  std::unordered_set<std::string> _snrExcludedWfs;
};

DEFINE_SMARTPOINTER(BatchLoader);

class BatchLoader : public Loader
{
public:
  BatchLoader(const std::string &recordStream);

  virtual ~BatchLoader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev);

  void load();

protected:
  bool _dataLoaded;
  IO::RecordStreamPtr _rs;
  std::unordered_multimap<std::string,
                          std::pair<const Core::TimeWindow, TimeWindowBuffer>>
      _streamMap;
  std::unordered_map<std::string, GenericRecordCPtr> _waveforms;
};

} // namespace Waveform
} // namespace HDD
} // namespace Seiscomp

#endif
