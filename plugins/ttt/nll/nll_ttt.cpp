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

#define SEISCOMP_COMPONENT NLLGrid
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

#include "hdd/ttt/nllgrid.h"

using namespace std;
using namespace Seiscomp;
namespace fs = boost::filesystem;

namespace {

/**
 * NLLGrid
 *
 * A class to compute seismic travel times for an homogeneous (constant)
 * velocity
 */
class NLLGrid : public TravelTimeTableInterface
{
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
  std::unordered_set<std::string> _validPphases;
  std::unordered_set<std::string> _validSphases;
  std::unique_ptr<HDD::TTT::NLLGrid> _grids;
};

double computeDistance(double lat1,
                       double lon1,
                       double lat2,
                       double lon2,
                       double *azimuth     = nullptr,
                       double *backAzimuth = nullptr)
{
  double dist;
  Math::Geo::delazi(lat1, lon1, lat2, lon2, &dist, azimuth, backAzimuth);
  return Math::Geo::deg2km(dist);
}

bool NLLGrid::setModel(const string &model)
{

  // load global configuration
  auto app = Seiscomp::System::Application::Instance();
  const Config::Config *cfg;
  Config::Config tmp;

  if (app)
  {
    cfg = &app->configuration();
  }
  else
  {
    if (!Environment::Instance()->initConfig(&tmp, ""))
    {
      return false;
    }
    else
    {
      cfg = &tmp;
    }
  }

  string base = "ttt.nllgrid." + model + ".";

  std::string gridPath;
  std::string gridModel;
  bool swapBytes        = false;
  unsigned maxOpenFiles = 512;

  try
  {
    std::vector<std::string> toks = cfg->getStrings(base + "P-Phases");
    _validPphases =
        std::unordered_set<std::string>(std::make_move_iterator(toks.begin()),
                                        std::make_move_iterator(toks.end()));
  }
  catch (...)
  {
    _validPphases = {"P", "Pg", "Pn"};
  }

  try
  {
    std::vector<std::string> toks = cfg->getStrings(base + "S-Phases");
    _validPphases =
        std::unordered_set<std::string>(std::make_move_iterator(toks.begin()),
                                        std::make_move_iterator(toks.end()));
  }
  catch (...)
  {
    _validSphases = {"S", "Sg", "Sn"};
  }

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

  try
  {
    swapBytes = cfg->getBool(base + "swapBytes");
  }
  catch (...)
  {}

  try
  {
    maxOpenFiles = cfg->getInt(base + "maxOpenFiles");
  }
  catch (...)
  {}

  _grids.reset(
      new HDD::TTT::NLLGrid(gridPath, gridModel, swapBytes, maxOpenFiles));

  _model = model;
  return true;
}

const string &NLLGrid::model() const { return _model; }

TravelTimeList *NLLGrid::compute(double lat1,
                                 double lon1,
                                 double dep1,
                                 double lat2,
                                 double lon2,
                                 double alt2,
                                 int ellc)
{
  TravelTimeList *ttlist = new TravelTimeList;
  ttlist->delta          = computeDistance(lat1, lon1, lat2, lon2);
  ttlist->depth          = dep1;
  try
  {
    ttlist->push_back(compute("P", lat1, lon1, dep1, lat2, lon2, alt2, ellc));
  }
  catch (const NoPhaseError &e)
  {}
  try
  {
    ttlist->push_back(compute("S", lat1, lon1, dep1, lat2, lon2, alt2, ellc));
  }
  catch (const NoPhaseError &e)
  {}
  ttlist->sortByTime();
  return ttlist;
}

TravelTime NLLGrid::compute(const char *phase,
                            double lat1,
                            double lon1,
                            double dep1,
                            double lat2,
                            double lon2,
                            double alt2,
                            int ellc)
{
  std::string phaseType;
  if (_validPphases.count(phase) > 0)
  {
    phaseType = "P";
  }
  else if (_validSphases.count(phase) > 0)
  {
    phaseType = "S";
  }
  else
  {
    throw NoPhaseError();
  }

  try
  {
    double travelTime, azimuth, takeOffAngle, velocity;
    _grids->compute(lat1, lon1, dep1, lat2, lon2, alt2, phaseType, travelTime,
                    azimuth, takeOffAngle, velocity);
    double dtdd = std::cos(takeOffAngle) // [sec/deg]
                  / Math::Geo::km2deg(velocity);
    double dtdh = std::sin(takeOffAngle) / velocity; // [sec/km]
    return TravelTime(phase, travelTime, dtdd, dtdh, 0, takeOffAngle);
  }
  catch (HDD::Exception &e)
  {
    SEISCOMP_WARNING("%s", e.what());
    throw NoPhaseError();
  }
}

double NLLGrid::computeTime(const char *phase,
                            double lat1,
                            double lon1,
                            double dep1,
                            double lat2,
                            double lon2,
                            double alt2,
                            int ellc)
{
  std::string phaseType;
  if (_validPphases.count(phase) > 0)
  {
    phaseType = "P";
  }
  else if (_validSphases.count(phase) > 0)
  {
    phaseType = "S";
  }
  else
  {
    throw NoPhaseError();
  }

  try
  {
    return _grids->compute(lat1, lon1, dep1, lat2, lon2, alt2, phaseType);
  }
  catch (HDD::Exception &e)
  {
    SEISCOMP_WARNING("%s", e.what());
    throw NoPhaseError();
  }
}

TravelTime NLLGrid::computeFirst(double lat1,
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
    "NonLinLoc Grid Travel Time Tables", "Luca Scarabello, ETH Zurich", 1, 0, 0)
REGISTER_TRAVELTIMETABLE(NLLGrid, "tttnll");

} // namespace
