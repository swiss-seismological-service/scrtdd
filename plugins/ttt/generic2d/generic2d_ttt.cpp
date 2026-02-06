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

#define SEISCOMP_COMPONENT Generic2D
#include <seiscomp/logging/log.h>

#include <seiscomp/core/plugin.h>
#include <seiscomp/core/strings.h>
#include <seiscomp/datamodel/config.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/seismology/ttt.h>
#include <seiscomp/system/application.h>
#include <seiscomp/system/environment.h>

#include <boost/filesystem.hpp>
#include <string>

#include "hdd/log.h"
#include "hdd/ttt/generic2d.h"

using namespace std;
using namespace Seiscomp;
namespace fs = boost::filesystem;

namespace {

/**
 * Generic2D
 * A class to read 2D seismic travel times from csv files
 */
class Generic2D : public TravelTimeTableInterface
{
public:
  Generic2D() { initializeLibraryOnce(); }

public:
  bool setModel(const std::string &model) override;
  const std::string &model() const override;

  TravelTimeList *compute(double lat1,
                          double lon1,
                          double dep1,
                          double lat2,
                          double lon2,
                          double alt2 = 0.,
                          int ellc    = 1) override;

  TravelTime compute(const char *phase,
                     double lat1,
                     double lon1,
                     double dep1,
                     double lat2,
                     double lon2,
                     double alt2 = 0.,
                     int ellc    = 1) override;

  TravelTime computeFirst(double lat1,
                          double lon1,
                          double dep1,
                          double lat2,
                          double lon2,
                          double alt2 = 0.,
                          int ellc    = 1) override;

  double computeTime(const char *phase,
                     double lat1,
                     double lon1,
                     double dep1,
                     double lat2,
                     double lon2,
                     double elev2 = 0.,
                     int ellc     = 1) override;

private:
  std::string _model;
  HDD::TTT::Generic2D _gen2d;

private:
  static bool initializeLibraryOnce()
  {
    static bool initialized = []() {
      // init HDD Logger
      auto hddLogger = [](HDD::Logger::Level level, const std::string &msg) {
        switch (level)
        {
        case HDD::Logger::Level::debug: SEISCOMP_DEBUG_S(msg); break;
        case HDD::Logger::Level::info: SEISCOMP_INFO_S(msg); break;
        case HDD::Logger::Level::warning: SEISCOMP_WARNING_S(msg); break;
        case HDD::Logger::Level::error: SEISCOMP_ERROR_S(msg); break;
        case HDD::Logger::Level::none: break;
        }
      };
      HDD::Logger::setLogger(hddLogger);
      HDD::Logger::setLevel(HDD::Logger::Level::debug);
      return true;
    }();
    return initialized;
  }
};

bool Generic2D::setModel(const string &model)
{
  _model = "";
  _gen2d = std::move(HDD::TTT::Generic2D());

  auto app = Seiscomp::System::Application::Instance();
  const Config::Config *cfg;
  Config::Config tmp;

  if (app) // app specific configuration
  {
    cfg = &app->configuration();
  }
  else // load global configuration
  {
    if (!Environment::Instance()->initConfig(&tmp, ""))
    {
      return false;
    }
    cfg = &tmp;
  }

  string base = "ttt.Generic2D." + model + ".";

  std::string gridPath;
  std::string gridModel;

  try
  {
    fs::path p = Environment::Instance()->absolutePath(
        cfg->getString(base + "tablePath"));
    gridPath  = p.parent_path().string();
    gridModel = p.filename().string();
  }
  catch (...)
  {
    return false;
  }

  _gen2d = std::move(HDD::TTT::Generic2D(gridPath, gridModel));

  _model = model;
  return true;
}

const string &Generic2D::model() const { return _model; }

TravelTimeList *Generic2D::compute(double lat1,
                                   double lon1,
                                   double dep1,
                                   double lat2,
                                   double lon2,
                                   double alt2,
                                   int ellc)
{
  double dist;
#if SC_API_VERSION < SC_API_VERSION_CHECK(16, 0, 0)
  double az, baz;
  Math::Geo::delazi(lat1, lon1, lat2, lon2, &dist, &az, &baz);
#else
  Math::Geo::delazi(lat1, lon1, lat2, lon2, &dist);
#endif

  TravelTimeList *ttlist = new TravelTimeList;
  ttlist->delta          = dist;
  ttlist->depth          = dep1;

  for (const string &phase : _gen2d.availablePhases())
  {
    try
    {
      ttlist->push_back(
          compute(phase.c_str(), lat1, lon1, dep1, lat2, lon2, alt2, ellc));
    }
    catch (const NoPhaseError &e)
    {}
  }
  ttlist->sortByTime();
  return ttlist;
}

TravelTime Generic2D::compute(const char *phase,
                              double lat1,
                              double lon1,
                              double dep1,
                              double lat2,
                              double lon2,
                              double alt2,
                              int ellc)
{
  try
  {
    TravelTime tt(phase, 0, 0, 0, 0, 0);
    double takeOffAzi;
    _gen2d.compute(lat1, lon1, dep1, lat2, lon2, alt2, phase, tt.time,
                   takeOffAzi, tt.takeoff, tt.dtdd, tt.dtdh);
#if SC_API_VERSION >= SC_API_VERSION_CHECK(16, 0, 0)
    tt.azi = takeOffAzi;
#endif
    return tt;
  }
  catch (HDD::TTT::Generic2D::OutOfRange &e)
  {
    throw std::out_of_range(e.what());
  }
  catch (HDD::Exception &e)
  {
    throw NoPhaseError();
  }
}

double Generic2D::computeTime(const char *phase,
                              double lat1,
                              double lon1,
                              double dep1,
                              double lat2,
                              double lon2,
                              double alt2,
                              int ellc)
{
  try
  {
    return _gen2d.compute(lat1, lon1, dep1, lat2, lon2, alt2, phase);
  }
  catch (HDD::TTT::Generic2D::OutOfRange &e)
  {
    throw std::out_of_range(e.what());
  }
  catch (HDD::Exception &e)
  {
    throw NoPhaseError();
  }
}

TravelTime Generic2D::computeFirst(double lat1,
                                   double lon1,
                                   double dep1,
                                   double lat2,
                                   double lon2,
                                   double alt2,
                                   int ellc)
{
  std::unique_ptr<TravelTimeList> ttlist(
      compute(lat1, lon1, dep1, lat2, lon2, alt2, ellc));
  if (!ttlist || ttlist->empty())
  {
    throw NoPhaseError();
  }
  return ttlist->front();
}

ADD_SC_PLUGIN(
    "Generic2D Travel Time Tables", "Luca Scarabello, ETH Zurich", 1, 0, 0)
REGISTER_TRAVELTIMETABLE(Generic2D, "Generic2D");

} // namespace
