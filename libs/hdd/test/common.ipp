
#include "catalog.h"
#include "utils.h"

#include <vector>

#include <boost/filesystem.hpp>

namespace {

struct TTTParams
{
  std::string type;
  std::string model;
};

const std::vector<TTTParams> _tttList = {
    {"LOCSAT", "iasp91"},
    {"libtau", "iasp91"},
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
     "./data/nll/iasp91_3D_simple/model/iasp91.PHASE.mod;"
     "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.time;"
     "./data/nll/iasp91_3D_simple/time/iasp91.PHASE.STATION.angle"},
    {"NonLinLoc",
     "./data/nll/iasp91_3D_sdc/model/iasp91.PHASE.mod;"
     "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.time;"
     "./data/nll/iasp91_3D_sdc/time/iasp91.PHASE.STATION.angle"}
};

// Some nll grids are too big to be added to git repository and need to
// be generated (see: hdd/test/data/nll/generate.sh) before running the
// tests. We check here what is available
std::vector<TTTParams> getTTTList()
{
  std::vector<TTTParams> ttts;
  for (const auto &prms : _tttList)
  {
    if (prms.type == "NonLinLoc")
    {
      static const std::regex regex(";", std::regex::optimize);
      std::vector<std::string> tokens(Seiscomp::HDD::splitString(prms.model, regex));
      const boost::filesystem::path velGridPath   = tokens.at(0);
      const boost::filesystem::path timeGridPath  = tokens.at(1);
      const boost::filesystem::path angleGridPath = tokens.at(2);
      if (!boost::filesystem::is_empty(velGridPath.parent_path())  ||
          !boost::filesystem::is_empty(timeGridPath.parent_path()) ||
          !boost::filesystem::is_empty(angleGridPath.parent_path()))
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

const std::vector<TTTParams> tttList = getTTTList();

// Those station parameters must be consistent with nonlinloc
// grids control files
// For tests with locsat it doesn't matter
std::vector<Seiscomp::HDD::Catalog::Station> stationList = {
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

