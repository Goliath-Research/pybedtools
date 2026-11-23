#ifndef LINEFILEUTILITIES_H
#define LINEFILEUTILITIES_H

#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace std;

// templated function to convert objects to strings
template <typename T>
inline
std::string ToString(const T & value) {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

// Split str on any character in delimiter, skipping empty tokens and stopping
// at the first NUL (the same tokens strtok yields). Works on str in place:
// lines come from input files and may be arbitrarily long, so they must never
// be copied into a stack buffer sized from their length.
template <typename Append>
inline
void TokenizeInPlace(const string &str, const string &delimiter, Append append)
{
    const string delims(delimiter.c_str());
    string::size_type end = str.find('\0');
    if (end == string::npos) end = str.size();
    string::size_type start = str.find_first_not_of(delims, 0);
    while (start != string::npos && start < end) {
        string::size_type stop = str.find_first_of(delims, start);
        if (stop == string::npos || stop > end) stop = end;
        append(str, start, stop - start);
        start = str.find_first_not_of(delims, stop);
    }
}

struct AppendString {
    vector<string> &elems;
    explicit AppendString(vector<string> &e) : elems(e) {}
    void operator()(const string &s, string::size_type pos, string::size_type len) const {
        elems.push_back(s.substr(pos, len));
    }
};

struct AppendInt {
    vector<int> &elems;
    explicit AppendInt(vector<int> &e) : elems(e) {}
    void operator()(const string &s, string::size_type pos, string::size_type len) const {
        elems.push_back(atoi(s.substr(pos, len).c_str()));
    }
};

// tokenize into a list of strings.
inline
void Tokenize(const string &str, vector<string> &elems, const string &delimiter = "\t") 
{
    TokenizeInPlace(str, delimiter, AppendString(elems));
}

// tokenize into a list of integers
inline
void Tokenize(const string &str, vector<int> &elems, const string &delimiter = "\t") {
    TokenizeInPlace(str, delimiter, AppendInt(elems));
}

#endif /* LINEFILEUTILITIES_H */

