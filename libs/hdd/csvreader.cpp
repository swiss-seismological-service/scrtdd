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
#include <fstream>
#include <iterator>

#include "csvreader.h"

using namespace std;

namespace HDD {
namespace CSV {

enum class CSVState
{
  UnquotedField,
  QuotedField,
  QuotedQuote
};

vector<string> readRow(const string &row)
{
  CSVState state = CSVState::UnquotedField;
  vector<string> fields{""};
  size_t i = 0; // index of the current field
  for (char c : row)
  {
    switch (state)
    {
    case CSVState::UnquotedField:
      switch (c)
      {
      case ',': // end of field
        fields.push_back("");
        i++;
        break;
      case '"': state = CSVState::QuotedField; break;
      default: fields[i].push_back(c); break;
      }
      break;
    case CSVState::QuotedField:
      switch (c)
      {
      case '"': state = CSVState::QuotedQuote; break;
      default: fields[i].push_back(c); break;
      }
      break;
    case CSVState::QuotedQuote:
      switch (c)
      {
      case ',': // , after closing quote
        fields.push_back("");
        i++;
        state = CSVState::UnquotedField;
        break;
      case '"': // "" -> "
        fields[i].push_back('"');
        state = CSVState::QuotedField;
        break;
      default: // end of quote
        state = CSVState::UnquotedField;
        break;
      }
      break;
    }
  }
  return fields;
}

// Read CSV file, Excel dialect.
vector<vector<string>> read(istream &in)
{
  vector<vector<string>> table;
  string row;
  while (!in.eof())
  {
    getline(in, row);
    if (in.bad() || in.fail())
    {
      break;
    }
    // handle windows end-of-line
    if (!row.empty() && *row.rbegin() == '\r')
    {
      row.pop_back();
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

vector<unordered_map<string, string>>
format(const vector<string> &header,
       const vector<vector<string>>::const_iterator &begin,
       const vector<vector<string>>::const_iterator &end)
{
  vector<unordered_map<string, string>> rows;
  for (auto it = begin; it != end; it++)
  {
    const vector<string> columns = *it;
    unordered_map<string, string> row;
    for (size_t i = 0; i < header.size(); ++i)
    {
      row[header[i]] = i < columns.size() ? columns[i] : "";
    }
    rows.push_back(row);
  }
  return rows;
}

vector<unordered_map<string, string>> readWithHeader(istream &in)
{
  vector<vector<string>> rows = read(in);
  return format(rows[0], rows.begin() + 1, rows.end());
}

vector<unordered_map<string, string>>
readWithHeader(istream &in, const vector<string> &header)
{
  vector<vector<string>> rows = read(in);
  return format(header, rows.begin(), rows.end());
}

vector<unordered_map<string, string>> readWithHeader(const string &filename)
{
  ifstream csvfile;
  csvfile.exceptions(std::ios::failbit | std::ios::badbit);
  csvfile.open(filename);
  csvfile.exceptions(std::ios::goodbit);
  return readWithHeader(csvfile);
}

vector<unordered_map<string, string>>
readWithHeader(const string &filename, const vector<string> &header)
{
  ifstream csvfile;
  csvfile.exceptions(std::ios::failbit | std::ios::badbit);
  csvfile.open(filename);
  csvfile.exceptions(std::ios::goodbit);
  return readWithHeader(csvfile, header);
}

} // namespace CSV
} // namespace HDD
