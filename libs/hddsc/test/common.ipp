
#include "hdd/catalog.h"
#include "hdd/nllttt.h"
#include "hdd/cvttt.h"
#include "hdd/utils.h"
#include "scttt.h"

#include <vector>

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace {

struct TTTParams
{
  std::string type;
  std::string model;
};

// Some nll grids are too big to be added to git repository and need to
// be generated (see: hdd/test/data/nll/generate.sh) before running the
// tests. We check here what is available
std::vector<TTTParams> __getTTTList__()
{
  static const std::vector<TTTParams> tttList = {
      {"LOCSAT", "iasp91"},
      {"libtau", "iasp91"},
      {"homogeneous", "iasp91"},
      {"ConstVel", "5.8;3.36"},
      {"NonLinLoc",
       "./data/nll/iasp91_2D_simple/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_2D_simple/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_2D_simple/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
        "./data/nll/iasp91_2D_sdc/model/iasp91.PHASE.mod;"
        "./data/nll/iasp91_2D_sdc/time/iasp91.PHASE.STATION.time;"
        "./data/nll/iasp91_2D_sdc/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
       "./data/nll/iasp91_2D_global/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_2D_global/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_2D_global/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
       "./data/nll/iasp91_2D_azimuthal_equidist/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_2D_azimuthal_equidist/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_2D_azimuthal_equidist/time/iasp91.PHASE.STATION.angle"}, 
      {"NonLinLoc",
       "./data/nll/iasp91_2D_merc/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_2D_merc/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_2D_merc/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
       "./data/nll/iasp91_2D_lambert/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_2D_lambert/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_2D_lambert/time/iasp91.PHASE.STATION.angle"}, 
      {"NonLinLoc",
       "./data/nll/iasp91_3D_simple/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
       "./data/nll/iasp91_3D_sdc/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
       "./data/nll/iasp91_3D_azimuthal_equidist/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_3D_azimuthal_equidist/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_3D_azimuthal_equidist/time/iasp91.PHASE.STATION.angle"}, 
      {"NonLinLoc",
       "./data/nll/iasp91_3D_merc/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_3D_merc/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_3D_merc/time/iasp91.PHASE.STATION.angle"},
      {"NonLinLoc",
       "./data/nll/iasp91_3D_lambert/model/iasp91.PHASE.mod;"
       "./data/nll/iasp91_3D_lambert/time/iasp91.PHASE.STATION.time;"
       "./data/nll/iasp91_3D_lambert/time/iasp91.PHASE.STATION.angle"}};

  std::vector<TTTParams> ttts;
  for (const auto &prms : tttList)
  {
    if (prms.type == "NonLinLoc")
    {
      static const std::regex regex(";", std::regex::optimize);
      std::vector<std::string> tokens(HDD::splitString(prms.model, regex));

      const fs::path velGridPath   = tokens.at(0);
      const fs::path timeGridPath  = tokens.at(1);
      const fs::path angleGridPath = tokens.at(2);

      if (!HDD::directoryEmpty(velGridPath.parent_path().string()) ||
          !HDD::directoryEmpty(timeGridPath.parent_path().string()) ||
          !HDD::directoryEmpty(angleGridPath.parent_path().string()))
      {
        ttts.push_back(prms);
      }
    }
    else
    {
      ttts.push_back(prms);
    }
  }
  return ttts;
}

const std::vector<TTTParams> tttList(__getTTTList__());

std::unique_ptr<HDD::TravelTimeTable> createTTT(const TTTParams &prms)
{
  if (prms.type == "NonLinLoc")
  {
    static const std::regex regex(";", std::regex::optimize);
    std::vector<std::string> tokens(HDD::splitString(prms.model, regex));

    const fs::path velGridPath   = tokens.at(0);
    const fs::path timeGridPath  = tokens.at(1);
    const fs::path angleGridPath = tokens.at(2);

    if (!HDD::directoryEmpty(velGridPath.parent_path().string()) ||
        !HDD::directoryEmpty(timeGridPath.parent_path().string()) ||
        !HDD::directoryEmpty(angleGridPath.parent_path().string()))
    {
      return std::unique_ptr<HDD::TravelTimeTable>(
          new HDD::NLL::TravelTimeTable(velGridPath.string(),
                                        timeGridPath.string(),
                                        angleGridPath.string(), false));
    }
  }
  else if (prms.type == "ConstVel")
  {
    static const std::regex regex(";", std::regex::optimize);
    std::vector<std::string> tokens(HDD::splitString(prms.model, regex));

    double pVel = std::stod(tokens.at(0));
    double sVel = std::stod(tokens.at(1));
    return std::unique_ptr<HDD::TravelTimeTable>(
        new HDD::ConstantVelocity(pVel, sVel));
  }
  else
  {
    return std::unique_ptr<HDD::TravelTimeTable>(
        new HDD::SCAdapter::TravelTimeTable(prms.type, prms.model));
  }
  return nullptr;
}

// Those station parameters must be consistent with nonlinloc
// grids control files
// For tests with locsat it doesn't matter
const std::vector<HDD::Catalog::Station> stationList = {
    {"NET.ST01A", 47.1, 8.6, 250, "NET", "ST01A", ""},
    {"NET.ST02A", 47.1, 8.4, 295, "NET", "ST02A", ""},
    {"NET.ST03A", 46.9, 8.4, 301, "NET", "ST03A", ""},
    {"NET.ST04A", 46.9, 8.6, 395, "NET", "ST04A", ""},
    {"NET.ST01B", 47.0, 8.7, 212, "NET", "ST01B", ""},
    {"NET.ST02B", 47.0, 8.3, 346, "NET", "ST02B", ""},
    {"NET.ST03B", 47.2, 8.5, 351, "NET", "ST03B", ""},
    {"NET.ST04B", 46.8, 8.5, 268, "NET", "ST04B", ""},
};

} // namespace
