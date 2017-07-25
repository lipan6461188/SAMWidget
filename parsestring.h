#ifndef PARSESTRING_H
#define PARSESTRING_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// split a string with a char
vector<string> split(const string &String, const char &c);

// a higher efficient method to split s a string
void split(const string &String, const char &c, vector<string> &Vec);

// trimming a string
void trim(string &String, const char &c);


#endif /* ParseString_hpp */
