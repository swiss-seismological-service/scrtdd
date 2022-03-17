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
#include "trace.h"
#include "utils.h"

#include <unordered_map>
#include <unordered_set>

namespace HDD {
namespace Waveform {

class Loader
{
public:
  Loader()          = default;
  virtual ~Loader() = default;

  Loader(const Loader &other) = delete;
  Loader &operator=(const Loader &other) = delete;

  // shared_ptr allows internal caching (unique_ptr would not) and
  // that's also why Trace is const
  virtual std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                           const Catalog::Phase &ph) = 0;
};

class BasicLoader : public Loader
{

public:
  BasicLoader(const std::string &recordStream) : _recordStreamURL(recordStream)
  {}
  virtual ~BasicLoader() = default;

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph) override;

  // ugly, but we need to give the user some feedbacks
  unsigned _counters_wf_no_avail   = 0;
  unsigned _counters_wf_downloaded = 0;

private:
  const std::string _recordStreamURL;
};

class BatchLoader : public Loader
{
public:
  BatchLoader(const std::string &recordStream)
      : _recordStreamURL(recordStream), _dataLoaded(false)
  {}

  virtual ~BatchLoader() = default;

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph) override;

  void request(const TimeWindow &tw, const Catalog::Phase &ph);

  size_t numRequests() { return _requests.size(); }

  void load();

  // ugly, but we need to give the user some feedbacks
  unsigned _counters_wf_no_avail   = 0;
  unsigned _counters_wf_downloaded = 0;

private:
  const std::string _recordStreamURL;
  bool _dataLoaded;
  std::unordered_multimap<std::string, const TimeWindow> _requests;
  std::unordered_map<std::string, std::shared_ptr<const Trace>> _traces;
};

class ExtraLenLoader : public Loader
{
public:
  ExtraLenLoader(const std::shared_ptr<Loader> &auxLdr, double traceMinLen)
      : ExtraLenLoader(auxLdr, traceMinLen / 2, traceMinLen / 2)
  {}

  ExtraLenLoader(const std::shared_ptr<Loader> &auxLdr,
                 double beforePickLen,
                 double afterPickLen)
      : _auxLdr(auxLdr), _beforePickLen(beforePickLen),
        _afterPickLen(afterPickLen)
  {}

  virtual ~ExtraLenLoader() = default;

  void setAuxLoader(const std::shared_ptr<Loader> &auxLdr) { _auxLdr = auxLdr; }
  std::shared_ptr<Loader> getAuxLoader() const { return _auxLdr; }

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph) override;

  TimeWindow traceTimeWindowToLoad(const TimeWindow &neededTW,
                                   const UTCTime &pickTime) const;

private:
  std::shared_ptr<Loader> _auxLdr;
  double _beforePickLen; // secs
  double _afterPickLen;  // secs
};

class DiskCachedLoader : public Loader
{
public:
  DiskCachedLoader(const std::shared_ptr<Loader> &auxLdr,
                   const std::string &cacheDir)
      : _auxLdr(auxLdr), _cacheDir(cacheDir)
  {}

  virtual ~DiskCachedLoader() = default;

  void setAuxLoader(const std::shared_ptr<Loader> &auxLdr) { _auxLdr = auxLdr; }
  std::shared_ptr<Loader> getAuxLoader() const { return _auxLdr; }

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph) override;

  bool isCached(const TimeWindow &tw,
                const Catalog::Phase &ph,
                const Catalog::Event &ev) const;

  // ugly, but we need to give the user some feedbacks
  unsigned _counters_wf_cached = 0;

private:
  std::unique_ptr<Trace> getFromCache(const TimeWindow &tw,
                                      const std::string &networkCode,
                                      const std::string &stationCode,
                                      const std::string &locationCode,
                                      const std::string &channelCode);

  void storeInCache(const TimeWindow &tw,
                    const std::string &networkCode,
                    const std::string &stationCode,
                    const std::string &locationCode,
                    const std::string &channelCode,
                    const Trace &trace);

  std::string waveformPath(const std::string &cacheDir,
                           const TimeWindow &tw,
                           const std::string &networkCode,
                           const std::string &stationCode,
                           const std::string &locationCode,
                           const std::string &channelCode) const;

  std::shared_ptr<Loader> _auxLdr;
  std::string _cacheDir;
};

class Processor
{
public:
  enum class Transform
  {
    NONE,
    L2,
    TRANSVERSAL,
    RADIAL
  };

  Processor()          = default;
  virtual ~Processor() = default;

  Processor(const Processor &other) = delete;
  Processor &operator=(const Processor &other) = delete;

  // shared_ptr allows internal caching (unique_ptr would not) and
  // that's also why Trace is const
  virtual std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                           const Catalog::Phase &ph,
                                           const Catalog::Event &ev,
                                           const Catalog::Station &sta,
                                           const std::string &filterStr,
                                           double resampleFreq,
                                           Transform trans) = 0;
};

class BasicProcessor : public Processor
{
public:
  BasicProcessor(const std::shared_ptr<Loader> &auxLdr) : _auxLdr(auxLdr) {}

  virtual ~BasicProcessor() = default;

  void setAuxLoader(const std::shared_ptr<Loader> &auxLdr) { _auxLdr = auxLdr; }
  std::shared_ptr<Loader> getAuxLoader() const { return _auxLdr; }

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const Catalog::Event &ev,
                                   const Catalog::Station &sta,
                                   const std::string &filterStr,
                                   double resampleFreq,
                                   Transform trans) override;

private:
  std::shared_ptr<Trace> process(const Trace &trace,
                                 const std::string &filterStr,
                                 double resampleFreq) const;
  std::shared_ptr<Loader> _auxLdr;
};

class MemCachedProc : public Processor
{
public:
  MemCachedProc(const std::shared_ptr<Processor> &auxPrc) : _auxPrc(auxPrc) {}

  virtual ~MemCachedProc() = default;

  void setAuxProcessor(const std::shared_ptr<Processor> &auxPrc)
  {
    _auxPrc = auxPrc;
  }
  std::shared_ptr<Processor> getAuxProcessor() const { return _auxPrc; }

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const Catalog::Event &ev,
                                   const Catalog::Station &sta,
                                   const std::string &filterStr,
                                   double resampleFreq,
                                   Transform trans) override;

  bool isCached(const TimeWindow &tw,
                const Catalog::Phase &ph,
                const Catalog::Event &ev) const;

private:
  std::shared_ptr<const Trace> getFromCache(const std::string &wfId);

  void storeInCache(const std::string &wfId,
                    const std::shared_ptr<const Trace> &trace);

  std::shared_ptr<Processor> _auxPrc;
  std::unordered_map<std::string, std::shared_ptr<const Trace>> _traces;
  std::unordered_set<std::string> _unloadables;
};

class SnrFilterPrc : public Processor
{

public:
  SnrFilterPrc(const std::shared_ptr<Processor> &auxPrc,
               double minSnr,
               double noiseStart,
               double noiseEnd,
               double signalStart,
               double signalEnd)
      : _auxPrc(auxPrc), _snr{minSnr, noiseStart, noiseEnd, signalStart,
                              signalEnd}
  {}

  virtual ~SnrFilterPrc() = default;

  void setAuxProcessor(const std::shared_ptr<Processor> &auxPrc)
  {
    _auxPrc = auxPrc;
  }
  std::shared_ptr<Processor> getAuxProcessor() const { return _auxPrc; }

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const Catalog::Event &ev,
                                   const Catalog::Station &sta,
                                   const std::string &filterStr,
                                   double resampleFreq,
                                   Transform trans) override;

  bool enabled() { return _enabled; }
  void setEnabled(bool enabled) { _enabled = enabled; }

  TimeWindow snrTimeWindow(const UTCTime &pickTime) const;

  bool goodSnr(const Trace &trace, const UTCTime &pickTime) const;

  unsigned _counters_wf_snr_low = 0;

private:
  std::shared_ptr<Processor> _auxPrc;
  struct
  {
    double minSnr;
    double noiseStart;  // secs relative to pick time
    double noiseEnd;    // secs relative to pick time
    double signalStart; // secs relative to pick time
    double signalEnd;   // secs relative to pick time
  } _snr;
  bool _enabled = true;
};

struct ThreeComponents
{
  enum Component
  {
    Vertical         = 0, /* usually Z */
    FirstHorizontal  = 1, /* usually N */
    SecondHorizontal = 2  /* usually E */
  };
  std::string names[3];
  double gain[3];
  double dip[3];
  double azimuth[3];
};

inline std::string getBandAndInstrumentCodes(const std::string &channelCode)
{
  if (channelCode.size() >= 2) return channelCode.substr(0, 2);
  return "";
}

inline std::string getOrientationCode(const std::string &channelCode)
{
  if (channelCode.size() == 3) return channelCode.substr(2, 3);
  return "";
}

void writeTrace(const Trace &trace, const std::string &file);
std::unique_ptr<Trace> readTrace(const std::string &file);

void resample(Trace &trace, double new_sf);

void filter(Trace &trace,
            bool demeaning,
            const std::string &filterStr,
            double resampleFreq);

double computeSnr(const Trace &tr,
                  const UTCTime &pickTime,
                  double noiseOffsetStart,
                  double noiseOffsetEnd,
                  double signalOffsetStart,
                  double signalOffsetEnd);

std::unique_ptr<Trace> transformL2(const TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const ThreeComponents &tc,
                                   const Trace trH1,
                                   const Trace trH2);

std::unique_ptr<Trace> transformRT(const TimeWindow &tw,
                                   const Catalog::Phase &ph,
                                   const Catalog::Event &ev,
                                   const Catalog::Station &sta,
                                   const ThreeComponents &tc,
                                   const Trace trV,
                                   const Trace trH1,
                                   const Trace trH2,
                                   Processor::Transform trans);

} // namespace Waveform
} // namespace HDD

#endif
