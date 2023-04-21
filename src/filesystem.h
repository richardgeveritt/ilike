//
//  filesystem.h
//  ilike_cpp
//
//  Created by Richard Everitt on 23/11/2022.
//  Code from https://stackoverflow.com/questions/20358455/cross-platform-way-to-make-a-directory-including-subfolders

#ifndef FILESYSTEM_H
#define FILESYSTEM_H

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#if defined(_WIN32)
#include <direct.h>
#endif

inline bool directory_exists(const std::string &directory_name)
{
  struct stat info;
  
  int statRC = stat( directory_name.c_str(), &info );
  if( statRC != 0 )
  {
    if (errno == ENOENT)  { return 0; } // something along the path does not exist
    if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
    return false;
  }
  
  return ( info.st_mode & S_IFDIR ) ? true : false;
}

inline void make_directory(const std::string &sPath)
{
   mode_t nMode = 0733; // UNIX style permissions
   int nError = 0;
   #if defined(_WIN32)
   nError = _mkdir(sPath.c_str()); // can be used on Windows
   #else
   nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
   #endif
   if (nError != 0) {
     stop("Error creating directory.");
   }
}

#endif // #ifndef INCLUDE_STD_FILESYSTEM_EXPERIMENTAL
