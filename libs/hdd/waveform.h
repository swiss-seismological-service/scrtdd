/***************************************************************************
 * MIT License                                                             *
 *                                                                         *
 * Copyright (C) by ETHZ/SED                                               *
 *                                                                         *
 * Permission is hereby granted, free of charge, to any person obtaining a *
 * copy of this software and associated documentation files (the           *
 * “Software”), to deal in the Software without restriction, including     *
 * without limitation the rights to use, copy, modify, merge, publish,     *
 * distribute, sublicense, and/or sell copies of the Software, and to      *
 * permit persons to whom the Software is furnished to do so, subject to   *
 * the following conditions:                                               *
 *                                                                         *
 * The above copyright notice and this permission notice shall be          *
 * included in all copies or substantial portions of the Software.         *
 *                                                                         *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,         *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      *
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  *
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    *
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    *
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __HDD_WAVEFORM_H__
#define __HDD_WAVEFORM_H__

#include "catalog.h"
#include "trace.h"
#include "utils.h"

#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace HDD {
namespace Waveform {

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

class Proxy
{
public:
  Proxy()          = default;
  virtual ~Proxy() = default;

  Proxy(const Proxy &other)            = delete;
  Proxy &operator=(const Proxy &other) = delete;

  // load one trace at the requested time window for a specific station
  // throw std::exception is the trace cannot be loaded
  virtual std::unique_ptr<Trace> loadTrace(const TimeWindow &tw,
                                           const std::string &networkCode,
                                           const std::string &stationCode,
                                           const std::string &locationCode,
                                           const std::string &channelCode) = 0;

  // load multiple traces: this is used in single-event relocation and when
  // preloading all the catalog traces. The reason of this method is to avoid
  // multiple re-connections to the data source (e.g. fdsn, seedlink) when
  // loading multiple traces. If there is no overhead in the initial connection,
  // this method can be implemented by calling loadTrace multiple times. For
  // each succefully loaded trace the 'onTraceLoaded' callback is called,
  // otherwise 'onTraceFailed' is called
  virtual void loadTraces(
      const std::unordered_multimap<std::string, const TimeWindow> &request,
      const std::function<void(const std::string &,
                               const TimeWindow &,
                               std::unique_ptr<Trace>)> &onTraceLoaded,
      const std::function<void(const std::string &,
                               const TimeWindow &,
                               const std::string &)> &onTraceFailed) = 0;

  // Return orientation information for a station. This is used to
  // perform rotation on the RT components or to know which components
  // are horizontal or vertical.
  // If this information is not knwon then throw std::exception and the
  // code will still be able to work as long as the user doesn't configure
  // R, T, H, or V components in HDD::Config::xcorr::components
  virtual void getComponentsInfo(const Catalog::Phase &ph,
                                 ThreeComponents &components) = 0;

  // Run a bandpass filter on the trace. The syntax of the 'filterStr' is
  // implementation specific. The 'filterStr' is configured in
  // HDD::Config::wfFilter::filterStr and passed to this function
  // when the filtering is required
  // throw HDD::Exception when an error happens
  virtual void filter(Trace &trace, const std::string &filterStr) = 0;

  // Write and Read traces. Used to create the waveform cache and to dump
  // waveforms for debugging. Choosing a standard format like miniseed allows
  // the user to inspect the traces with an external tool
  // throw std::exception when an error happens
  virtual void writeTrace(const Trace &trace, const std::string &file) = 0;
  virtual std::unique_ptr<Trace> readTrace(const std::string &file)    = 0;
};

// This class can be used when no-crosscorrelation is required
class NoWaveformProxy : public Proxy
{
public:
  std::unique_ptr<HDD::Trace> loadTrace(const HDD::TimeWindow &tw,
                                        const std::string &networkCode,
                                        const std::string &stationCode,
                                        const std::string &locationCode,
                                        const std::string &channelCode) override
  {
    throw HDD::Exception("No waveform available");
  }

  void loadTraces(
      const std::unordered_multimap<std::string, const HDD::TimeWindow>
          &request,
      const std::function<void(const std::string &,
                               const HDD::TimeWindow &,
                               std::unique_ptr<HDD::Trace>)> &onTraceLoaded,
      const std::function<void(const std::string &,
                               const HDD::TimeWindow &,
                               const std::string &)> &onTraceFailed) override
  {}

  void getComponentsInfo(const HDD::Catalog::Phase &ph,
                         HDD::Waveform::ThreeComponents &components) override
  {
    throw HDD::Exception("No waveform available");
  }

  void filter(HDD::Trace &trace, const std::string &filterStr) override
  {
    throw HDD::Exception("Not implemented");
  }

  void writeTrace(const HDD::Trace &trace, const std::string &file) override
  {
    throw HDD::Exception("Not implemented");
  }

  std::unique_ptr<HDD::Trace> readTrace(const std::string &file) override
  {
    throw HDD::Exception("Not implemented");
  }
};

class Loader
{
public:
  Loader()          = default;
  virtual ~Loader() = default;

  Loader(const Loader &other)            = delete;
  Loader &operator=(const Loader &other) = delete;

  // shared_ptr allows internal caching (unique_ptr would not) and
  // that's also why Trace is const
  virtual std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                           const Catalog::Phase &ph) = 0;
};

class BasicLoader : public Loader
{

public:
  BasicLoader(const std::shared_ptr<Proxy> &wf) : _wf(wf) {}
  virtual ~BasicLoader() = default;

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph) override;

  // ugly, but we need to give the user some feedbacks
  unsigned _counters_wf_no_avail   = 0;
  unsigned _counters_wf_downloaded = 0;

private:
  std::shared_ptr<Proxy> _wf;
};

class BatchLoader : public Loader
{
public:
  BatchLoader(const std::shared_ptr<Proxy> &wf) : _wf(wf), _dataLoaded(false) {}

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
  std::shared_ptr<Proxy> _wf;
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
  DiskCachedLoader(const std::shared_ptr<Proxy> &wf,
                   const std::shared_ptr<Loader> &auxLdr,
                   const std::string &cacheDir)
      : _wf(wf), _auxLdr(auxLdr), _cacheDir(cacheDir)
  {}

  virtual ~DiskCachedLoader() = default;

  void setAuxLoader(const std::shared_ptr<Loader> &auxLdr) { _auxLdr = auxLdr; }
  std::shared_ptr<Loader> getAuxLoader() const { return _auxLdr; }

  std::shared_ptr<const Trace> get(const TimeWindow &tw,
                                   const Catalog::Phase &ph) override;

  bool isCached(const TimeWindow &tw,
                const Catalog::Phase &ph,
                const Catalog::Event &ev) const;

  void writeTrace(const Trace &trace, const std::string &file);
  std::unique_ptr<Trace> readTrace(const std::string &file);

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

  std::shared_ptr<Proxy> _wf;
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

  Processor(const Processor &other)            = delete;
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
  BasicProcessor(const std::shared_ptr<Proxy> &wf,
                 const std::shared_ptr<Loader> &auxLdr,
                 double extraTraceLen)
      : _wf(wf), _auxLdr(auxLdr), _extraTraceLen(extraTraceLen)
  {}

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

  void filter(Trace &trace,
              bool demeaning,
              const std::string &filterStr,
              double resampleFreq) const;

private:
  std::shared_ptr<Trace> loadAndProcess(const TimeWindow &tw,
                                        Catalog::Phase ph,
                                        const std::string &channelCode,
                                        const std::string &filterStr,
                                        double resampleFreq) const;
  std::shared_ptr<Proxy> _wf;
  std::shared_ptr<Loader> _auxLdr;
  double _extraTraceLen;
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

inline std::string getBandAndInstrumentCodes(const std::string &channelCode)
{
  if (channelCode.size() >= 2) return channelCode.substr(0, 2);
  return "";
}

inline std::string getOrientationCode(const std::string &channelCode)
{
  if (channelCode.size() == 3) return channelCode.substr(2, 1);
  return "";
}

void resample(Trace &trace, double new_sf);

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
