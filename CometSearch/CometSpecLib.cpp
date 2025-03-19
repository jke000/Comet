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
#include <string.h>


Mutex g_vSpecLibMutex;


CometSpecLib::CometSpecLib()
{
}


CometSpecLib::~CometSpecLib()
{
}


// SpecLib will be a vector of structs. Structs contain spectra and IDs (peptide, scan #, etc.)
bool CometSpecLib::LoadSpecLib(string strSpecLibFile)
{
   FILE *fp;
   size_t tTmp;
   char szBuf[SIZE_BUF];

   if (g_bSpecLibRead)
      return true;

   if ((fp = fopen(g_staticParams.speclibInfo.strSpecLibFile.c_str(), "r")) == NULL)
   {
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error, spectral library file cannot be read: '%s'.\n", g_staticParams.speclibInfo.strSpecLibFile.c_str());
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }
   else
      fclose(fp);

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

   printf("OK speclib extension %s\n", strExtension.c_str() );

   if (strExtension == ".db")
   {
      if (!ReadSpecLibSqlite(strSpecLibFile))
         return false;
   }
   else if (strExtension == ".raw")
   {
      if (!ReadSpecLibRaw(strSpecLibFile))
         return false;
   }
   else if (strExtension == ".msp")
   {
      if (!ReadSpecLibMSP(strSpecLibFile))
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

   for (auto it = g_vSpecLib.begin(); it != g_vSpecLib.end(); ++it)
   {
      printf("OK.  %d, %s, %d peaks\n", (*it).iLibEntry, (*it).strName.c_str(), (*it).iNumPeaks);
      
      for (int i = 0; i < (*it).iNumPeaks; ++i)
      {
         printf("\t%0.2lf\t%0.2lf\n", (*it).vSpecLibPeaks.at(i).first, (*it).vSpecLibPeaks.at(i).second);
         if (i == 4)
            break;
      }
   }

   return true;
}


bool CometSpecLib::ReadSpecLibSqlite(string strSpecLibFile)
{

   printf(" Error - sqlite/.db files as spectral libraries are not supported yet.\n");
   exit(1);

   sqlite3* db;
   sqlite3_stmt* stmt;
   const char* sql = "SELECT * FROM SpectrumTable";

   // Open the database
   if (sqlite3_open(strSpecLibFile.c_str(), &db) != SQLITE_OK)
   {
//    std::cerr << "Cannot open sqlite database: " << sqlite3_errmsg(db) << std::endl;
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error - cannot open sqlite database file '%s': %s.\n", strSpecLibFile.c_str(), sqlite3_errmsg(db));
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


bool CometSpecLib::ReadSpecLibRaw(string strSpecLibFile)
{
   printf(" Error - raw files as spectral libraries are not supported yet.\n");
   exit(1);

   MSReader mstReader;

   vector<MSSpectrumType> msLevel;

   mstReader.setFilter(msLevel);

   if (g_staticParams.options.iSpecLibMSLevel == 1)
      msLevel.push_back(MS1);
   else if (g_staticParams.options.iSpecLibMSLevel == 2)
      msLevel.push_back(MS2);
   else if (g_staticParams.options.iSpecLibMSLevel == 3)
      msLevel.push_back(MS3);
   else
   {
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error, MS level not set for the spectral library input.\n");
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

/*
   Spectrum mstSpectrum;           // For holding spectrum.

   // Load all input spectra.
   while (true)
   {
      // Loads in MSMS spectrum data.
      if (_bFirstScan)
      {
         PreloadIons(mstReader, mstSpectrum, false, 0);  // Use 0 as scan num here in last argument instead of iFirstScan; must
         _bFirstScan = false;                            // be MS/MS scan else data not read by MSToolkit so safer to start at 0.
      }                                                  // Not ideal as could be reading non-relevant scans but it's fast enough.
      else
      {
         PreloadIons(mstReader, mstSpectrum, true);
      }

      if (iFileLastScan == -1)
      {
         iFileLastScan = mstReader.getLastScan();
         if (iLastScan == 0 && iFileLastScan != -1)
            iLastScan = iFileLastScan;
      }

      if ((iFileLastScan != -1) && (iFileLastScan < iFirstScan))
      {
         _bDoneProcessingAllSpectra = true;
         break;
      }

      iScanNumber = mstSpectrum.getScanNumber();

      if (iLastScan > 0 && iScanNumber > iLastScan)
         break;

      if (g_staticParams.bSkipToStartScan && iScanNumber < iFirstScan)
      {
         g_staticParams.bSkipToStartScan = false;

         PreloadIons(mstReader, mstSpectrum, false, iFirstScan);
         iScanNumber = mstSpectrum.getScanNumber();

         // iScanNumber will equal 0 if iFirstScan is not the right scan level
         // So need to keep reading the next scan until we get a non-zero scan number
         while (iScanNumber == 0 && iFirstScan < iFileLastScan)
         {
            iFirstScan++;
            PreloadIons(mstReader, mstSpectrum, false, iFirstScan);
            iScanNumber = mstSpectrum.getScanNumber();
         }
      }
      if (iScanNumber != 0)
      {
         if (iLastScan > 0 && iScanNumber > iLastScan)
         {
            _bDoneProcessingAllSpectra = true;
            break;
         }
         if (iFirstScan != 0 && iLastScan > 0 && !(iFirstScan <= iScanNumber && iScanNumber <= iLastScan))
            continue;
         if (iFirstScan != 0 && iLastScan == 0 && iScanNumber < iFirstScan)
            continue;

         g_staticParams.bSkipToStartScan = false;
         iTmpCount = iScanNumber;

         // To run a search, all that's needed is MH+ and Z. So need to generate
         // all combinations of these for each spectrum, whether there's a known
         // Z for each precursor or if Comet has to guess the 1+ or 2+/3+ charges.

         for (int i = 0 ; i < mstSpectrum.sizeMZ(); ++i)  // walk through all precursor m/z's; usually just one
         {
            double dMZ = 0.0;              // m/z to use for analysis
            vector<int> vChargeStates;

            if (mstSpectrum.sizeMZ() == mstSpectrum.sizeZ())
            {
               iSpectrumCharge = mstSpectrum.atZ(i).z;
            }
            else if (mstSpectrum.sizeMZ() == 1 && mstSpectrum.sizeMZ() < mstSpectrum.sizeZ())
            {
               // example from ms2 file with one precursor and multiple Z lines?
               // will need to include all spectrum charges below
               iSpectrumCharge = mstSpectrum.atZ(i).z;
            }
            else if (mstSpectrum.sizeMZ() > mstSpectrum.sizeZ())
            {
               // need to ignore any spectrum charge as don't know which correspond charge to which precursor
               iSpectrumCharge = 0;
            }
            else
            {
               // don't know what condition/spectrum type leads here
               iSpectrumCharge = 0;
               printf(" Warning, scan %d has %d precursors and %d precursor charges\n", iScanNumber, mstSpectrum.sizeMZ(), mstSpectrum.sizeZ());
            }

            // Thermo's monoisotopic m/z determine can fail sometimes. Assume that when
            // the mono m/z value is less than selection window, it is wrong and use the
            // selection m/z as the precursor m/z. This should
            // be invoked when searching Thermo raw files and mzML converted from those.
            // Only applied when single precursor present.
            dMZ = mstSpectrum.getMonoMZ(i);

            if (g_staticParams.options.bCorrectMass && mstSpectrum.sizeMZ() == 1)
            {
               double dSelectionLower = mstSpectrum.getSelWindowLower();
               double dSelectedMZ = mstSpectrum.getMZ(i);

               if (dMZ > 0.1 && dSelectionLower > 0.1 && dMZ+0.1 < dSelectionLower)
                  dMZ = dSelectedMZ;
            }

                  for (int i=0; i<mstSpectrum.size(); ++i)
                  {
                     dSumTotal += mstSpectrum.at(i).intensity;

                     if (mstSpectrum.at(i).mz < mstSpectrum.getMZ())
                        dSumBelow += mstSpectrum.at(i).intensity;
                  }


*/

}


// MSP format:
//
// Name: AAAGELQEDSGLMALAK/2_0_30eV
// MW: 1675.8440
// Comment: Single Pep=Tryptic Mods=0 Fullname=K.AAAGELQEDSGLMALAK.L/2 Charge=2 Parent=837.9220 Mz_diff=1.0ppm  HCD=30.00% Scan=81146 Origfile="2018...
// Num peaks: 201
// 101.071	2366.7	"IQA/0.6ppm"
// 102.0552	4117.3	"IEA/2.4ppm"
// 105.0657	730.7	"?"
// 107.4682	612.6	"?"
bool CometSpecLib::ReadSpecLibMSP(string strSpecLibFile)
{
   FILE *fp;

   printf("OK in ReadSpecLibMSP\n");

   if ( (fp=fopen(strSpecLibFile.c_str(), "r")) == NULL)
   {
      char szErrorMsg[SIZE_ERROR];
      sprintf(szErrorMsg, "Error, MSP spectral library file cannot be read: '%s'.\n", strSpecLibFile.c_str());
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   char szBuf[SIZE_BUF];
   char szTmp[SIZE_BUF];
   int iWhichLibEntry = 0;

   fgets(szBuf, SIZE_BUF, fp);
   while (!feof(fp))
   {
      szBuf[SIZE_BUF - 1] = '\0'; // terminate string in case line was longer than SIZE_BUF


printf("OK szBuf %s\n", szBuf);
      if (!strncmp(szBuf, "Name:", 5))
      {
         iWhichLibEntry++;

printf("OK in Name string\n");

         if (szBuf[strlen(szBuf) - 1] != '\n') // really long Name: line, parse to newline char
         {
            char cChar;
            cChar = szBuf[strlen(szBuf)-1];
            while (cChar != '\n')
            {
               cChar = getc(fp);
            }
         }
         while (szBuf[strlen(szBuf) - 1] == '\r' || szBuf[strlen(szBuf) - 1] == '\n')
            szBuf[strlen(szBuf) - 1] = '\0';

         sscanf(szBuf + 6, "%s", szTmp);

         struct SpecLibStruct pTmp;
         pTmp.strName = szTmp;
         pTmp.iLibEntry = iWhichLibEntry;
         pTmp.iCharge = 0;
         pTmp.dMW = 0;

         while (fgets(szBuf, SIZE_BUF, fp))
         {
            if (!strncmp(szBuf, "Name:", 5))
            {
               break;
            }
            else if (!strncmp(szBuf, "MW:", 3))
            {
               // FIX:  in my example, the MW: field encodes (parent m/z * charge)
               sscanf(szBuf+4, "%lf", &(pTmp.dMW));
            }
            else if (!strncmp(szBuf, "Comment:", 8))
            {
               char *pStr;

               if ((pStr = strstr(szBuf, " Charge=")) != NULL)
                  sscanf(pStr + 9, "%d", &(pTmp.iCharge));
            }
            else if (!strncmp(szBuf, "Num peaks:", 10))
            {
               int iNumPeaks = 0;
               sscanf(szBuf + 11, "%d", &iNumPeaks);

               pTmp.iNumPeaks = iNumPeaks;

               for (int i = 0; i < iNumPeaks; ++i)
               {
                  double dMass;
                  double dInten;

                  fgets(szBuf, SIZE_BUF, fp);
                  sscanf(szBuf, "%lf %lf %*s", &dMass, &dInten);

                  if (dMass > 0.0 && dMass <1e6 && dInten > 0.0)  // some sanity check on parsed mass
                     pTmp.vSpecLibPeaks.push_back(make_pair(dMass, dInten));
               }

               break;
            }
         }

         pTmp.dMW -= pTmp.iCharge * PROTON_MASS;  // make neutral mass


         // FIX:  do something with the peak list depending on score/processing
         g_vSpecLib.push_back(pTmp);
printf("OK pushing g_vSpecLib size %d, NumPeaks %d\n", g_vSpecLib.size(), pTmp.iNumPeaks );

         SetSpecLibPrecursorIndex(pTmp.dMW, pTmp.iCharge, g_vSpecLib.size() - 1);

printf("OK after SetSpecLibPrecursorIndex\n");
      }
      else
         fgets(szBuf, SIZE_BUF, fp);

   }

   fclose(fp);

   return true;
}


double CometSpecLib::ScoreSpecLib(Query *pQuery,
                                  unsigned int iWhichSpecLib)
{
   double dScore = 0.0;
   int x, y, bin;

printf("OK in ScoreSpecLib, query %d, iWhichSpecLib %d\n", pQuery->_spectrumInfoInternal.iScanNumber, iWhichSpecLib);

   int iMax = pQuery->_spectrumInfoInternal.iArraySize/SPARSE_MATRIX_SIZE;

   for (auto it = g_vSpecLib.at(iWhichSpecLib).vSpecLibPeaks.begin(); it != g_vSpecLib.at(iWhichSpecLib).vSpecLibPeaks.end() ; ++it)
   {
      x = BIN(it->first);

      if (!(bin <= 0 || x>iMax || pQuery->ppfSparseFastXcorrData[x] == NULL))
      {
         y = bin - (x * SPARSE_MATRIX_SIZE);
         dScore += pQuery->ppfSparseFastXcorrData[x][y];
      }
   }

   dScore = std::round(dScore * 0.005 * 1000.0) / 1000.0;  // round to 3 decimal points

printf("OK in ScoreSpecLib, query %d, iWhichSpecLib %d, dScore %lf\n", pQuery->_spectrumInfoInternal.iScanNumber, iWhichSpecLib, dScore);

   return dScore;
}


// For each spec lib mass "bin", set g_vulSpecLibPrecursorIndex which is a vector of all
// SpecLib entries that are matched to that "bin". This allows a mass query to walk through
// and score against all entries in the vector.
void CometSpecLib::SetSpecLibPrecursorIndex(double dNeutralMass,
                                            int iCharge,
                                            size_t iWhichSpecLib)
{
   printf("OK in SetSpecLibPrecursorIndex\n");

   double dProtonatedMass = dNeutralMass + PROTON_MASS;

   double dToleranceLow = 0;
   double dToleranceHigh = 0;

   int iMaxBin = BINPREC(g_staticParams.options.dPeptideMassHigh);

   printf("OK2 in SetSpecLibPrecursorIndex, iWhichSpecLib %ld, g_vSpecLib.size() %ld\n", iWhichSpecLib, g_vSpecLib.size());
   int iPrecursorCharge = g_vSpecLib.at(iWhichSpecLib).iCharge;
   printf("OK3 in SetSpecLibPrecursorIndex\n");

   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      dToleranceLow  = g_staticParams.tolerances.dInputToleranceMinus;
      dToleranceHigh = g_staticParams.tolerances.dInputTolerancePlus;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         dToleranceLow  *= iPrecursorCharge;
         dToleranceHigh *= iPrecursorCharge;
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      dToleranceLow  = g_staticParams.tolerances.dInputToleranceMinus * 0.001;
      dToleranceHigh = g_staticParams.tolerances.dInputTolerancePlus  * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         dToleranceLow  *= iPrecursorCharge;
         dToleranceHigh *= iPrecursorCharge;
      }
   }

   // tolerances are fixed above except if ppm is specified
   else if (g_staticParams.tolerances.iMassToleranceUnits == 2) // ppm
   {
      dToleranceLow  = g_staticParams.tolerances.dInputToleranceMinus * dProtonatedMass / 1E6;
      dToleranceHigh = g_staticParams.tolerances.dInputTolerancePlus * dProtonatedMass / 1E6;
   }
   else
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg,  " Error - peptide_mass_units must be 0, 1 or 2. Value set is %d.\n", g_staticParams.tolerances.iMassToleranceUnits);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   // these are the range of neutral mass bins if any theoretical peptide falls into,
   // we want to add them to the fragment index
   double dMassLow = dProtonatedMass + dToleranceLow;
   double dMassHigh = dProtonatedMass + dToleranceHigh;
   int iStart = BINPREC(dMassLow);  // add dToleranceLow as it will be negative number
   int iEnd   = BINPREC(dMassHigh);

   if (iStart < 0)
      iStart = 0;   // real problems if we actually get here
   if (iEnd > iMaxBin)
      iEnd = iMaxBin;

printf("*** OK iStart %d, iEnd %d, dMassLow %lf, dMassHigh %lf\n", iStart, iEnd, dMassLow, dMassHigh);
   for (int x = iStart ; x <= iEnd; ++x)
      g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);

printf("OK2 iStart %d, iEnd %d, dMassLow %lf, dMassHigh %lf\n", iStart, iEnd, dMassLow, dMassHigh);
   // now go through each isotope offset
   if (g_staticParams.tolerances.iIsotopeError > 0)
   {
      if (g_staticParams.tolerances.iIsotopeError >= 1
            && g_staticParams.tolerances.iIsotopeError <= 6)
      {
         iStart = BINPREC(dMassLow - C13_DIFF * PROTON_MASS);     // do +1 offset
         iEnd   = BINPREC(dMassHigh - C13_DIFF * PROTON_MASS);
         if (iStart < 0)
            iStart = 0;
         if (iEnd > iMaxBin)
            iEnd = iMaxBin;
         for (int x = iStart ; x <= iEnd; ++x)
         {
printf("OK1 pushback %d\n", x);
            g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
         }

         if (g_staticParams.tolerances.iIsotopeError >= 2
               && g_staticParams.tolerances.iIsotopeError <= 6
               && g_staticParams.tolerances.iIsotopeError != 5)
         {
            iStart = BINPREC(dMassLow - 2.0 * C13_DIFF * PROTON_MASS);     // do +2 offset
            iEnd   = BINPREC(dMassHigh - 2.0 * C13_DIFF * PROTON_MASS);
            if (iStart < 0)
               iStart = 0;
            if (iEnd > iMaxBin)
               iEnd = iMaxBin;
            for (int x = iStart ; x <= iEnd; ++x)
            {
printf("OK2 pushback %d\n", x);
               g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
            }

            if (g_staticParams.tolerances.iIsotopeError >= 3
                  && g_staticParams.tolerances.iIsotopeError <= 6
                  && g_staticParams.tolerances.iIsotopeError != 5)
            {
               iStart = BINPREC(dMassLow - 3.0 * C13_DIFF * PROTON_MASS);     // do +3 offset
               iEnd   = BINPREC(dMassHigh - 3.0 * C13_DIFF * PROTON_MASS);
               if (iStart < 0)
                  iStart = 0;
               if (iEnd > iMaxBin)
                  iEnd = iMaxBin;
               for (int x = iStart ; x <= iEnd; ++x)
               {
printf("OK3 pushback %d\n", x);
                  g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
               }
            }
         }
      }

      if (g_staticParams.tolerances.iIsotopeError == 5
            || g_staticParams.tolerances.iIsotopeError == 6)
      {
         iStart = BINPREC(dMassLow + C13_DIFF * PROTON_MASS);      // do -1 offset
         iEnd   = BINPREC(dMassHigh + C13_DIFF * PROTON_MASS);
         if (iStart < 0)
            iStart = 0;
         if (iEnd > iMaxBin)
            iEnd = iMaxBin;
         for (int x = iStart ; x <= iEnd; ++x)
         {
printf("OK4 pushback %d\n", x);
            g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
         }

         if (g_staticParams.tolerances.iIsotopeError == 6)     // do -2 and -3 offsets
         {
            iStart = BINPREC(dMassLow + 2.0 * C13_DIFF * PROTON_MASS);
            iEnd   = BINPREC(dMassHigh + 2.0 * C13_DIFF * PROTON_MASS);
            if (iStart < 0)
               iStart = 0;
            if (iEnd > iMaxBin)
               iEnd = iMaxBin;
            for (int x = iStart ; x <= iEnd; ++x)
            {
printf("OK5 pushback %d\n", x);
               g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
            }

            iStart = BINPREC(dMassLow + 3.0 * C13_DIFF * PROTON_MASS);
            iEnd   = BINPREC(dMassHigh + 3.0 * C13_DIFF * PROTON_MASS);
            if (iStart < 0)
               iStart = 0;
            if (iEnd > iMaxBin)
               iEnd = iMaxBin;
            for (int x = iStart ; x <= iEnd; ++x)
            {
printf("OK6 pushback %d\n", x);
               g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
            }
         }
      }
      else if (g_staticParams.tolerances.iIsotopeError == 7)            // do -8, -4, +4, +8 offsets
      {
         iStart = BINPREC(dMassLow + 8.0 * C13_DIFF * PROTON_MASS);
         iEnd   = BINPREC(dMassHigh + 8.0 * C13_DIFF * PROTON_MASS);
         if (iStart < 0)
            iStart = 0;
         if (iEnd > iMaxBin)
            iEnd = iMaxBin;
         for (int x = iStart ; x <= iEnd; ++x)
         {
printf("OK7 pushback %d\n", x);
            g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
         }

         iStart = BINPREC(dMassLow + 4.0 * C13_DIFF * PROTON_MASS);
         iEnd   = BINPREC(dMassHigh + 4.0 * C13_DIFF * PROTON_MASS);
         if (iStart < 0)
            iStart = 0;
         if (iEnd > iMaxBin)
            iEnd = iMaxBin;
         for (int x = iStart ; x <= iEnd; ++x)
         {
printf("OK8 pushback %d\n", x);
            g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
         }

         iStart = BINPREC(dMassLow - 8.0 * C13_DIFF * PROTON_MASS);
         iEnd   = BINPREC(dMassHigh - 8.0 * C13_DIFF * PROTON_MASS);
         if (iStart < 0)
            iStart = 0;
         if (iEnd > iMaxBin)
            iEnd = iMaxBin;
         for (int x = iStart ; x <= iEnd; ++x)
         {
printf("OK9 pushback %d\n", x);
            g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
         }

         iStart = BINPREC(dMassLow - 4.0 * C13_DIFF * PROTON_MASS);
         iEnd   = BINPREC(dMassHigh - 4.0 * C13_DIFF * PROTON_MASS);
         if (iStart < 0)
            iStart = 0;
         if (iEnd > iMaxBin)
            iEnd = iMaxBin;
         for (int x = iStart ; x <= iEnd; ++x)
         {
printf("OK10 pushback %d\n", x);
            g_vulSpecLibPrecursorIndex.at(x).push_back(iWhichSpecLib);
         }
      }
   }
}


void CometSpecLib::StoreSpecLib(Query *it,
                                unsigned int iWhichSpecLib,
                                double dSpecLibScore)
{

}


// Function to decode BLOB data as an array of 8-byte floats
std::vector<double> CometSpecLib::decodeBlob(const void* blob, int size)
{
   std::vector<double> result;

   const double* data = static_cast<const double*>(blob);

   int numDoubles = size / sizeof(double);

   for (int i = 0; i < numDoubles; ++i)
      result.push_back(data[i]);

   return result;
}


// Function to print the contents of a vector of doubles
void CometSpecLib::printDoubleVector(const std::vector<double>& vec)
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
