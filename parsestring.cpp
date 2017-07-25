//
//  ParseString.cpp
//  calcRPKM
//
//  Created by Lee on 2017/7/22.
//  Copyright © 2017年 Lee. All rights reserved.
//

#include "parsestring.h"

// split a string with a char
vector<string> split(const string &String, const char &c)
{
    vector<string> Vec;
    size_t last = 0;
    size_t begin = String.find(c, 0);
    while(begin != string::npos)
    {
        Vec.emplace_back(String.substr(last, begin-last));
        last = begin + 1;
        begin = String.find(c, begin+1);
    }
    Vec.emplace_back(String.substr(last, begin-last));
    return Vec;
}

// a higher efficient method to split s a string
void split(const string &String, const char &c, vector<string> &Vec)
{
    size_t last = 0;
    Vec.clear();
    size_t begin = String.find(c, 0);
    while(begin != string::npos)
    {
        Vec.emplace_back(String.substr(last, begin-last));
        last = begin + 1;
        begin = String.find(c, begin+1);
    }
    Vec.emplace_back(String.substr(last, begin-last));
}

// trimming a string
void trim(string &String, const char &c)
{
    size_t start_point = 0;
    for(;start_point<String.size();start_point++)
        if(String.at(start_point) != c)
            break;

    size_t end_point = String.size() - 1;
    for(;end_point>=start_point;end_point--)
        if(String.at(end_point) != c)
            break;

    if(end_point < String.size())
        String.erase(end_point+1, String.size());

    if(start_point > 0)
        String.erase(0, start_point);
}

