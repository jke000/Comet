// Copyright 2023 Jimmy Eng
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include "Common.h"
#include "CometSpecLib.h"
#include "CometSearch.h"
#include "ThreadPool.h"
#include "CometStatus.h"
#include "CometMassSpecUtils.h"
#include <string>


Mutex g_vSpecLibMutex;


CometSpecLib::CometSpecLib()
{
}


CometSpecLib::~CometSpecLib()
{
}


// SpecLib will be a vector of structs. Structs contain spectra and IDs (peptide, scan #, etc.)
bool CometSpecLib::ReadSpecLib(string strSpecLibFile)
{
   FILE *fp;
   size_t tTmp;
   char szBuf[SIZE_BUF];

   if (g_bSpecLibRead)
      return true;

   // Transform to lower case for case insensitive file extension match
   std::string strLowerFileName = strSpecLibFile;
   std::transform(strLowerFileName.begin(), strLowerFileName.end(), strLowerFileName.begin(), ::tolower);

   // Find the position of the last dot in the filename
   size_t dotPos = strLowerFileName.rfind('.');
   if (dotPos == std::string::npos)
   {
      // No dot found, so no extension
      return false;
   }

   // Extract the file extension
   std::string strExtension = strLowerFileName.substr(dotPos);

   if (strExtension == ".db")
   {
      if (!ReadSpecLibSqlite(g_staticParams.speclibInfo.strSpecLib))
         return false;
   }
   else if (strExtension == ".raw")
   {
      if (!ReadSpecLibRaw(g_staticParams.speclibInfo.strSpecLib))
         return false;
   }
   else
   {
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error, expecting sqlite .db or Thermo .raw file for the spectral library.\n");
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   g_bSpecLibRead = true;

   return true;
}


bool CometSpecLib::ReadSpecLibSqlite(string strSpecLib)
{
   sqlite3* db;
   sqlite3_stmt* stmt;
   const char* sql = "SELECT * FROM SpectrumTable";

   // Open the database
   if (sqlite3_open(strSpecLib.c_str(), &db) != SQLITE_OK)
   {
//    std::cerr << "Cannot open sqlite database: " << sqlite3_errmsg(db) << std::endl;
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error - cannot open sqlite database file '%s': %s.\n", strSpecLib.c_str(), sqlite3_errmsg(db));
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   // Prepare the SQL statement
   if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK)
   {
//    std::cerr << "Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error - sqlite failed to prepare statment: %s.\n", sqlite3_errmsg(db));
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      sqlite3_close(db);
      return false;
   }

   // Execute the SQL statement and read the rows
   while (sqlite3_step(stmt) == SQLITE_ROW)
   {
      int numCols = sqlite3_column_count(stmt);

      for (int col = 0; col < numCols; ++col)
      {
         const char* colName = sqlite3_column_name(stmt, col);
         int colType = sqlite3_column_type(stmt, col);

         std::cout << colName << ": ";

         switch (colType)
         {
            case SQLITE_INTEGER:
               std::cout << sqlite3_column_int(stmt, col) << std::endl;
               break;
            case SQLITE_FLOAT:
               std::cout << sqlite3_column_double(stmt, col) << std::endl;
               break;
            case SQLITE_TEXT:
               std::cout << sqlite3_column_text(stmt, col) << std::endl;
               break;
            case SQLITE_BLOB:
            {
               const void* blobData = sqlite3_column_blob(stmt, col);
               int blobSize = sqlite3_column_bytes(stmt, col);
               std::vector<double> decodedBlob = decodeBlob(blobData, blobSize);
               printDoubleVector(decodedBlob);
               break;
            }
            case SQLITE_NULL:
               std::cout << "NULL" << std::endl;
               break;
            default:
               std::cout << "Unknown data type" << std::endl;
               break;
         }
      }
      std::cout << "--------------------------------------" << std::endl;
   }

   // Clean up
   sqlite3_finalize(stmt);
   sqlite3_close(db);

   return 0;
}


bool CometSpecLib::ReadSpecLibRaw(string strSpecLib)
{
}


bool CometSpecLib::SearchSpecLib(int iWhichQuery,
                                 ThreadPool *tp)
{
   int i;

/*
   for (int iWhichLib = 0 ; iWhichLib < SpecLib.size() ; ++iWhichLib)
   {
      if (WithinSpecLibMassTolerance(iWhichQuery, iWhichLib))
      {
         SpecLibScore(iWhichQuery, iWhichLib);
      }
   }

*/
   return true;
}

// Function to decode BLOB data as an array of 8-byte floats
std::vector<double> decodeBlob(const void* blob, int size)
{
   std::vector<double> result;

   const double* data = static_cast<const double*>(blob);

   int numDoubles = size / sizeof(double);

   for (int i = 0; i < numDoubles; ++i)
      result.push_back(data[i]);

   return result;
}


// Function to print the contents of a vector of doubles
void printDoubleVector(const std::vector<double>& vec)
{
    std::cout << "[";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        std::cout << vec[i];

        if (i < vec.size() - 1)
        {
            std::cout << ", ";
        }

        if (i == 9)
        {
           printf(" ...");
           break;
        }
    }
    std::cout << "]" << std::endl;
}

