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

#ifndef __HDD_CSVREADER_H__
#define __HDD_CSVREADER_H__

#include <istream>
#include <string>
#include <unordered_map>
#include <vector>

namespace HDD {
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


 * Optionally, a header (i.e. the first line) can be present. If no header is
 * present, you can pass the header format to the function.
 * If a line has more fields than header fields, the excessing fields will
 * be discarded.
 * If a line has less fields than header fields, the missing ones will be set
 * to an empty string.

 header1,header2,header3,header4
 val1,val2,val3,val4
 val1,val2,val3,val4
 val1,val2,val3,val4
 val1,val2,val3,val4

 *
 */

/*
 * Read file without header line.
 */
std::vector<std::vector<std::string>> read(std::istream &in);

std::vector<std::vector<std::string>> read(const std::string &filename);

/*
 * Read file with header: The first line must be the header.
 */
std::vector<std::unordered_map<std::string, std::string>>
readWithHeader(std::istream &in);

std::vector<std::unordered_map<std::string, std::string>>
readWithHeader(const std::string &filename);

/*
 * Read file without header. The header is passed to the function instead.
 */
std::vector<std::unordered_map<std::string, std::string>>
readWithHeader(std::istream &in, const std::vector<std::string> &header);

std::vector<std::unordered_map<std::string, std::string>>
readWithHeader(const std::string &filename,
               const std::vector<std::string> &header);

} // namespace CSV
} // namespace HDD

#endif
