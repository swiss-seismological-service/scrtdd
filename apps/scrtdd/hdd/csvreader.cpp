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

#include <string>
#include <vector>
#include <iterator>
#include <istream>
#include <fstream>
#include <iostream>

#include "csvreader.h"

using namespace std;

namespace Seiscomp {
namespace HDD {
namespace CSV {

enum class CSVState {
    UnquotedField,
    QuotedField,
    QuotedQuote
};

vector<string> readRow(const string &row)
{
    CSVState state = CSVState::UnquotedField;
    vector<string> fields {""};
    size_t i = 0; // index of the current field
    for (char c : row)
    {
        switch (state)
        {
            case CSVState::UnquotedField:
                switch (c)
                {
                    case ',': // end of field
                              fields.push_back(""); i++;
                              break;
                    case '"': state = CSVState::QuotedField;
                              break;
                    default:  fields[i].push_back(c);
                              break;
                }
                break;
            case CSVState::QuotedField:
                switch (c)
                {
                    case '"': state = CSVState::QuotedQuote;
                              break;
                    default:  fields[i].push_back(c);
                              break;
                }
                break;
            case CSVState::QuotedQuote:
                switch (c)
                {
                    case ',': // , after closing quote
                              fields.push_back(""); i++;
                              state = CSVState::UnquotedField;
                              break;
                    case '"': // "" -> "
                              fields[i].push_back('"');
                              state = CSVState::QuotedField;
                              break;
                    default:  // end of quote
                              state = CSVState::UnquotedField;
                              break;
                }
                break;
        }
    }
    return fields;
}

/// Read CSV file, Excel dialect. Accept "quoted fields ""with quotes"""
vector<vector<string>> read(istream &in)
{
    vector<vector<string>> table;
    string row;
    while (!in.eof())
    {
        getline(in, row);
        if (in.bad() || in.fail()) {
            break;
        }
        auto fields = readRow(row);
        table.push_back(fields);
    }
    return table;
}

vector<vector<string>> read(const string &filename)
{
    ifstream csvfile(filename);
    return read(csvfile);
}

vector< unordered_map<string,string> > format(const vector<string>& header,
                                    const vector<vector<string>>::const_iterator& begin,
                                    const vector<vector<string>>::const_iterator& end)
{
    vector< unordered_map<string,string> > rows;
    for (auto it = begin; it != end; it++)
    {
        const vector<string> columns = *it;
        unordered_map<string,string> row;
        for (size_t i = 0; i < header.size(); ++i)
        {
            row[ header[i] ] = i < columns.size() ? columns[i] : "";
        }
        rows.push_back(row);
    }
    return rows;
}

vector< unordered_map<string,string> > readWithHeader(istream &in)
{
    vector<vector<string>> rows = read(in);
    return format(rows[0], rows.begin()+1, rows.end());
}

vector< unordered_map<string,string> > readWithHeader(istream &in,
                                            const vector<string>& header)
{
    vector<vector<string>> rows = read(in);
    return format(header, rows.begin(), rows.end());
}

vector< unordered_map<string,string> > readWithHeader(const string &filename)
{
    ifstream csvfile;
    csvfile.exceptions(std::ios::failbit | std::ios::badbit);
    csvfile.open(filename);
    csvfile.exceptions(std::ios::goodbit);
    return readWithHeader(csvfile);
}

vector< unordered_map<string,string> > readWithHeader(const string &filename,
                                            const vector<string>& header)
{
    ifstream csvfile;
    csvfile.exceptions(std::ios::failbit | std::ios::badbit);
    csvfile.open(filename);
    csvfile.exceptions(std::ios::goodbit);
    return readWithHeader(csvfile, header);
}

}
}
}
