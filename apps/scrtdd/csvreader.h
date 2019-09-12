/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the "SED Public License for Seiscomp Contributions"          *
 *                                                                         *
 *   You should have received a copy of the "SED Public License for        *
 *   Seiscomp Contributions" with this. If not, you can find it at         *
 *   http://www.seismo.ethz.ch/static/seiscomp_contrib/license.txt         *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   "SED Public License for Seiscomp Contributions" for more details      *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/


#ifndef __RTDD_APPLICATIONS_CSVREADER_H__
#define __RTDD_APPLICATIONS_CSVREADER_H__

namespace CSV {

/*
 * Read CSV file, Excel dialect
 *
 * E.g.

 val1,val2,val3,val4
 val1,val2,val3,val4
 val1,val2,val3,val4
 val1,val2,val3,val4

 * Accept "quoted fields ""with quotes"""
 * Each line can have a different number of columns
 *
 * E.g.

 val1,val2,val3,val4,val5,val6
 val1,"val 2 quoted",val3
 val1,"val 2 quoted with ""inner"" quote",val3,val4,val5
 val1,val2,"val 3 quoted, with separator",val4


 * Optionally a header, the first line, can be present. If not present you
 * can pass the header format to the function.
 * Each line has to have the same number of fields of the header, additional
 * columns will be disregarded

 header1,header2,header3,header4
 val1,val2,val3,val4
 val1,val2,val3,val4
 val1,val2,val3,val4
 val1,val2,val3,val4

 *
 */


/*
 * Read file without a header line
 */
std::vector<std::vector<std::string> > read(std::istream &in);

std::vector<std::vector<std::string> > read(const std::string &filename);

/*
 * Read file with a header: the first line must be the header
 */
std::vector< std::map<std::string,std::string> > readWithHeader(std::istream &in);

std::vector< std::map<std::string,std::string> > readWithHeader(const std::string &filename);

/*
 * Read file without a header. The header is passed to the function instead.
 */
std::vector< std::map<std::string,std::string> > readWithHeader(std::istream &in, const std::vector<std::string>& header);

std::vector< std::map<std::string,std::string> > readWithHeader(const std::string &filename, const std::vector<std::string>& header);

}

#endif
