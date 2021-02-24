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

#ifndef __HDD_WAVEFORM_H__
#define __HDD_WAVEFORM_H__

#include "catalog.h"

#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/datamodel/utils.h>

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
                             const std::string &channelCode);

bool merge(GenericRecord &trace, const RecordSequence &seq);
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
  Loader(const std::string &recordStream) : Loader(recordStream, false, false)
  {}

  Loader(LoaderPtr auxLdr) : Loader(auxLdr, false, false) {}

  virtual ~Loader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev,
                                bool demeaning               = false,
                                const std::string &filterStr = "",
                                double resampleFreq          = 0);

  virtual bool isCached(const Core::TimeWindow &tw,
                        const Catalog::Phase &ph,
                        const Catalog::Event &ev)
  {
    return false;
  }

  // empty string disables debugging
  void setDebugDirectory(const std::string &directory = "")
  {
    _wfDebugDir = directory;
  }

  // counters
  unsigned _counters_wf_no_avail   = 0;
  unsigned _counters_wf_cached     = 0;
  unsigned _counters_wf_downloaded = 0;

protected:
  Loader(const std::string &recordStream, bool doCaching, bool cacheProcessed)
      : _recordStreamURL(recordStream), _doCaching(doCaching),
        _cacheProcessed(cacheProcessed)
  {}

  Loader(LoaderPtr auxLdr, bool doCaching, bool cacheProcessed)
      : _auxLdr(auxLdr), _doCaching(doCaching), _cacheProcessed(cacheProcessed)
  {
    if (_doCaching && !_cacheProcessed && _auxLdr->_cacheProcessed)
    {
      throw std::runtime_error(
          "Cannot set processed cache as auxiliary of an unprocessed cache");
    }
  }

  virtual GenericRecordCPtr getFromCache(const Core::TimeWindow &tw,
                                         const std::string &networkCode,
                                         const std::string &stationCode,
                                         const std::string &locationCode,
                                         const std::string &channelCode)
  {
    return nullptr;
  }

  virtual void storeInCache(const Core::TimeWindow &tw,
                            const std::string &networkCode,
                            const std::string &stationCode,
                            const std::string &locationCode,
                            const std::string &channelCode,
                            const GenericRecordCPtr &trace)
  {}

  GenericRecordPtr readAndProjectWaveform(const Core::TimeWindow &tw,
                                          const Catalog::Phase &ph,
                                          const Catalog::Event &ev);

  LoaderPtr _auxLdr;
  const std::string _recordStreamURL;
  const bool _doCaching;
  const bool _cacheProcessed;

  std::string _wfDebugDir;
};

DEFINE_SMARTPOINTER(DiskCachedLoader);

class DiskCachedLoader : public Loader
{
public:
  DiskCachedLoader(const std::string &recordStream,
                   bool cacheProcessed,
                   const std::string &cacheDir)
      : Loader(recordStream, true, cacheProcessed), _cacheDir(cacheDir)
  {}

  DiskCachedLoader(LoaderPtr auxLdr,
                   bool cacheProcessed,
                   const std::string &cacheDir)
      : Loader(auxLdr, true, cacheProcessed), _cacheDir(cacheDir)
  {}

  virtual ~DiskCachedLoader() {}

  virtual bool isCached(const Core::TimeWindow &tw,
                        const Catalog::Phase &ph,
                        const Catalog::Event &ev);

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

class MemCachedLoader : public Loader
{

public:
  MemCachedLoader(const std::string &recordStream, bool cacheProcessed)
      : Loader(recordStream, true, cacheProcessed)
  {}

  MemCachedLoader(LoaderPtr auxLdr, bool cacheProcessed)
      : Loader(auxLdr, true, cacheProcessed)
  {}

  virtual ~MemCachedLoader() {}

  virtual bool isCached(const Core::TimeWindow &tw,
                        const Catalog::Phase &ph,
                        const Catalog::Event &ev);

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

class ExtraLenLoader : public Loader
{
public:
  ExtraLenLoader(const std::string &recordStream, double traceMinLen)
      : Loader(recordStream), _traceMinLen(traceMinLen)
  {}

  ExtraLenLoader(LoaderPtr auxLdr, double traceMinLen)
      : Loader(auxLdr), _traceMinLen(traceMinLen)
  {}

  virtual ~ExtraLenLoader() {}

  virtual GenericRecordCPtr get(const Core::TimeWindow &tw,
                                const Catalog::Phase &ph,
                                const Catalog::Event &ev,
                                bool demeaning               = false,
                                const std::string &filterStr = "",
                                double resampleFreq          = 0);

  Core::TimeWindow traceTimeWindowToLoad(const Core::TimeWindow &neededTW,
                                         const Core::Time &pickTime) const;

protected:
  double _traceMinLen; // secs
};

DEFINE_SMARTPOINTER(SnrFilteredLoader);

class SnrFilteredLoader : public Loader
{

public:
  SnrFilteredLoader(const std::string &recordStream,
                    double minSnr,
                    double noiseStart,
                    double noiseEnd,
                    double signalStart,
                    double signalEnd)
      : Loader(recordStream), _snr{minSnr, noiseStart, noiseEnd, signalStart,
                                   signalEnd}
  {}

  SnrFilteredLoader(LoaderPtr auxLdr,
                    double minSnr,
                    double noiseStart,
                    double noiseEnd,
                    double signalStart,
                    double signalEnd)
      : Loader(auxLdr), _snr{minSnr, noiseStart, noiseEnd, signalStart,
                             signalEnd}
  {}

  virtual ~SnrFilteredLoader() {}

  GenericRecordCPtr get(const Core::TimeWindow &tw,
                        const Catalog::Phase &ph,
                        const Catalog::Event &ev,
                        bool demeaning               = false,
                        const std::string &filterStr = "",
                        double resampleFreq          = 0);

  Core::TimeWindow snrTimeWindow(const Core::Time &pickTime) const;

  bool goodSnr(const GenericRecordCPtr &trace,
               const Core::Time &pickTime) const;

  // counters
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

} // namespace Waveform
} // namespace HDD
} // namespace Seiscomp

#endif
