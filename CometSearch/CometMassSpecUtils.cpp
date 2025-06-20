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


///////////////////////////////////////////////////////////////////////////////
//  Implementations for generic mass spectrometry related utility functions.
///////////////////////////////////////////////////////////////////////////////

#include "Common.h"
#include "CometDataInternal.h"
#include "CometSearchManager.h"
#include "CometMassSpecUtils.h"
#include "CometSearch.h"
#include <inttypes.h>


double CometMassSpecUtils::GetFragmentIonMass(int iWhichIonSeries,
                                              int i,
                                              int ctCharge,
                                              double *pdAAforward,
                                              double *pdAAreverse)
{
   double dFragmentIonMass = 0.0;

   switch (iWhichIonSeries)
   {
      case ION_SERIES_B:
         dFragmentIonMass = pdAAforward[i];
         break;

      case ION_SERIES_Y:
         dFragmentIonMass = pdAAreverse[i];
         break;

      case ION_SERIES_A:
         dFragmentIonMass = pdAAforward[i] - g_staticParams.massUtility.dCO;
         break;

      case ION_SERIES_C:
         dFragmentIonMass = pdAAforward[i] + g_staticParams.massUtility.dNH3;
         break;

      case ION_SERIES_X:
         dFragmentIonMass = pdAAreverse[i] + g_staticParams.massUtility.dCOminusH2;
         break;

      case ION_SERIES_Z:
         dFragmentIonMass = pdAAreverse[i] - g_staticParams.massUtility.dNH2;
         break;

      case ION_SERIES_Z1:
         dFragmentIonMass = pdAAreverse[i] - g_staticParams.massUtility.dNH2 + Hydrogen_Mono;
         break;
   }

   return (dFragmentIonMass + (ctCharge - 1.0) * PROTON_MASS) / ctCharge;
}


void CometMassSpecUtils::AssignMass(double *pdAAMass,
                                    int bMonoMasses,
                                    double *dOH2)
{
   double H, O, C, N, S, Se;

   if (bMonoMasses) // monoisotopic masses
   {
      H = pdAAMass[(int)'h'] = Hydrogen_Mono; // hydrogen
      O = pdAAMass[(int)'o'] = Oxygen_Mono;  // oxygen
      C = pdAAMass[(int)'c'] = Carbon_Mono;   // carbon
      N = pdAAMass[(int)'n'] = Nitrogen_Mono;   // nitrogen
//    P = pdAAMass[(int)'p'] = 30.973762;    // phosphorus
      S = pdAAMass[(int)'s'] = 31.9720707;   // sulphur
      Se = pdAAMass[(int)'e'] = 79.9165196;  // selenium
   }
   else  // average masses
   {
      H = pdAAMass[(int)'h'] =  1.00794;
      O = pdAAMass[(int)'o'] = 15.9994;
      C = pdAAMass[(int)'c'] = 12.0107;
      N = pdAAMass[(int)'n'] = 14.0067;
//    P = pdAAMass[(int)'p'] = 30.973761;
      S = pdAAMass[(int)'s'] = 32.065;
      Se = pdAAMass[(int)'e'] = 78.96;
   }

   *dOH2 = H + H + O;

   pdAAMass[(int)'G'] = C*2  + H*3  + N   + O ;
   pdAAMass[(int)'A'] = C*3  + H*5  + N   + O ;
   pdAAMass[(int)'S'] = C*3  + H*5  + N   + O*2 ;
   pdAAMass[(int)'P'] = C*5  + H*7  + N   + O ;
   pdAAMass[(int)'V'] = C*5  + H*9  + N   + O ;
   pdAAMass[(int)'T'] = C*4  + H*7  + N   + O*2 ;
   pdAAMass[(int)'C'] = C*3  + H*5  + N   + O   + S ;
   pdAAMass[(int)'U'] = C*3  + H*5  + N   + O   + Se ;
   pdAAMass[(int)'L'] = C*6  + H*11 + N   + O ;
   pdAAMass[(int)'I'] = C*6  + H*11 + N   + O ;
   pdAAMass[(int)'N'] = C*4  + H*6  + N*2 + O*2 ;
   pdAAMass[(int)'D'] = C*4  + H*5  + N   + O*3 ;
   pdAAMass[(int)'Q'] = C*5  + H*8  + N*2 + O*2 ;
   pdAAMass[(int)'K'] = C*6  + H*12 + N*2 + O ;
   pdAAMass[(int)'E'] = C*5  + H*7  + N   + O*3 ;
   pdAAMass[(int)'M'] = C*5  + H*9  + N   + O   + S ;
   pdAAMass[(int)'H'] = C*6  + H*7  + N*3 + O ;
   pdAAMass[(int)'F'] = C*9  + H*9  + N   + O ;
   pdAAMass[(int)'R'] = C*6  + H*12 + N*4 + O ;
   pdAAMass[(int)'Y'] = C*9  + H*9  + N   + O*2 ;
   pdAAMass[(int)'W'] = C*11 + H*10 + N*2 + O ;

   pdAAMass[(int)'O'] = C*12  + H*19 + N*3 + O*2 ;

   pdAAMass[(int)'B'] = 0.0;
   pdAAMass[(int)'J'] = 0.0;
   pdAAMass[(int)'X'] = 0.0;
   pdAAMass[(int)'Z'] = 0.0;
}


// return a single protein name as a C char string
void CometMassSpecUtils::GetProteinName(FILE *fpfasta,
                                        comet_fileoffset_t lFilePosition,
                                        char *szProteinName)
{
   size_t tTmp;

   comet_fseek(fpfasta, lFilePosition, SEEK_SET);

   if (g_staticParams.iIndexDb)  //fragment ion or peptide index
   {
      long lSize;

      tTmp = fread(&lSize, sizeof(long), 1, fpfasta);
      vector<comet_fileoffset_t> vOffsets;
      for (long x = 0; x < lSize; ++x) // read file offsets
      {
         comet_fileoffset_t tmpoffset;
         tTmp = fread(&tmpoffset, sizeof(comet_fileoffset_t), 1, fpfasta);
         vOffsets.push_back(tmpoffset);
      }
      for (long x = 0; x < lSize; ++x) // read name from fasta
      {
         char szTmp[WIDTH_REFERENCE];
         comet_fseek(fpfasta, vOffsets.at(x), SEEK_SET);
         tTmp = fread(szTmp, sizeof(char)*WIDTH_REFERENCE, 1, fpfasta);
         sscanf(szTmp, "%511s", szProteinName);  // WIDTH_REFERENCE-1
         break;  //break here to only get first protein reference (out of lSize)
      }
   }
   else  //regular fasta database
   {
      fscanf(fpfasta, "%511s", szProteinName);  // WIDTH_REFERENCE-1
      szProteinName[511] = '\0';
   }
}


// return a single protein sequence as C++ string
void CometMassSpecUtils::GetProteinSequence(FILE *fpfasta,
                                            comet_fileoffset_t lFilePosition,
                                            string &strSeq)
{
   strSeq.clear();

   if (!g_staticParams.iIndexDb)  // works only for regular FASTA
   {
      int iTmpCh;

      comet_fseek(fpfasta, lFilePosition, SEEK_SET);

      // skip to end of description line
      while (((iTmpCh = getc(fpfasta)) != '\n') && (iTmpCh != '\r') && (iTmpCh != EOF));

      // load sequence
      while (((iTmpCh=getc(fpfasta)) != '>') && (iTmpCh != EOF))
      {
         if ('a'<=iTmpCh && iTmpCh<='z')
         {
            strSeq += iTmpCh - 32;  // convert toupper case so subtract 32 (i.e. 'A'-'a')
            g_staticParams.databaseInfo.uliTotAACount++;
         }
         else if ('A'<=iTmpCh && iTmpCh<='Z')
         {
            strSeq += iTmpCh;
            g_staticParams.databaseInfo.uliTotAACount++;
         }
         else if (iTmpCh == '*')  // stop codon
         {
            strSeq += iTmpCh;
         }
      }
   }
}


// return all matched protein names in a vector of strings
void CometMassSpecUtils::GetProteinNameString(FILE *fpfasta,
                                              int iWhichQuery,  // which search
                                              int iWhichResult, // which peptide within the search
                                              int iPrintTargetDecoy,    // 0 = target+decoys, 1=target only, 2=decoy only
                                              bool bReturnFullProteinString,   // 0 = return accession only, 1 = return full description line
                                              unsigned int *uiNumTotProteins,   // matched protein count
                                              vector<string>& vProteinTargets,  // the target protein names
                                              vector<string>& vProteinDecoys)   // the decoy protein names if applicable
{
   char szProteinName[512];

   int iLenDecoyPrefix = (int)strlen(g_staticParams.szDecoyPrefix);

   // FIX:  protein references is so convoluted with the restoration of peptide index.  This
   // seems to work now but definitely needs to be revisited to be cleaned up.
   // Look into lProteinFilePosition and lWhichProtein with Results struct.

   if (g_staticParams.iIndexDb == 2)  // peptide index
   {
      Results* pOutput;

      if (iPrintTargetDecoy != 2)
         pOutput = g_pvQuery.at(iWhichQuery)->_pResults;
      else
         pOutput = g_pvQuery.at(iWhichQuery)->_pDecoys;

      *uiNumTotProteins = (unsigned int)g_pvProteinsList.at(pOutput[iWhichResult].lProteinFilePosition).size();

      int iPrintDuplicateProteinCt = 0; // track # proteins, exit when at iMaxDuplicateProteins

      vector<string> vTmp;      // store decoy matches here to append at end

      for (auto it = pOutput[iWhichResult].pWhichProtein.begin(); it != pOutput[iWhichResult].pWhichProtein.end(); ++it)
      {
         for (auto it2 = g_pvProteinsList.at((*it).lWhichProtein).begin(); it2 != g_pvProteinsList.at((*it).lWhichProtein).end(); ++it2)
         {
            comet_fileoffset_t lEntry = (*it2);
            comet_fseek(fpfasta, lEntry, SEEK_SET);

            if (bReturnFullProteinString)
               fgets(szProteinName, 511, fpfasta);
            else
               fscanf(fpfasta, "%500s", szProteinName);

            szProteinName[500] = '\0';  // limit protein name strings to 500 chars

            // remove all terminating chars
            while ((szProteinName[strlen(szProteinName) - 1] == '\n') || (szProteinName[strlen(szProteinName) - 1] == '\r'))
               szProteinName[strlen(szProteinName) - 1] = '\0';

            if (!strncmp(szProteinName, g_staticParams.szDecoyPrefix, iLenDecoyPrefix))
               vTmp.push_back(szProteinName);
            else
               vProteinTargets.push_back(szProteinName);

            iPrintDuplicateProteinCt++;
            if (iPrintDuplicateProteinCt >= g_staticParams.options.iMaxDuplicateProteins)
               break;
         }
         if (iPrintDuplicateProteinCt >= g_staticParams.options.iMaxDuplicateProteins)
            break;
      }
      if (vTmp.size() > 0)      // append any decoy matches now; these would be decoys present in fasta
         vProteinTargets.insert(vProteinTargets.end(), vTmp.begin(), vTmp.end());


      vTmp.clear();
      for (auto it = pOutput[iWhichResult].pWhichDecoyProtein.begin(); it != pOutput[iWhichResult].pWhichDecoyProtein.end(); ++it)
      {
         for (auto it2 = g_pvProteinsList.at((*it).lWhichProtein).begin(); it2 != g_pvProteinsList.at((*it).lWhichProtein).end(); ++it2)
         {
            comet_fileoffset_t lEntry = (*it2);
            comet_fseek(fpfasta, lEntry, SEEK_SET);

            if (bReturnFullProteinString)
               fgets(szProteinName, 511, fpfasta);
            else
               fscanf(fpfasta, "%500s", szProteinName);

            szProteinName[500] = '\0';  // limit protein name strings to 500 chars

            // remove all terminating chars
            while ((szProteinName[strlen(szProteinName) - 1] == '\n') || (szProteinName[strlen(szProteinName) - 1] == '\r'))
               szProteinName[strlen(szProteinName) - 1] = '\0';

            if (!strncmp(szProteinName, g_staticParams.szDecoyPrefix, iLenDecoyPrefix))
               vTmp.push_back(szProteinName);
            else
               vProteinDecoys.push_back(szProteinName);

            iPrintDuplicateProteinCt++;
            if (iPrintDuplicateProteinCt >= g_staticParams.options.iMaxDuplicateProteins)
               break;
         }
         if (iPrintDuplicateProteinCt >= g_staticParams.options.iMaxDuplicateProteins)
            break;
      }
      if (vTmp.size() > 0)      // append any decoy matches now; these would be decoys present in fasta
         vProteinTargets.insert(vProteinTargets.end(), vTmp.begin(), vTmp.end());
   }
   else if (g_staticParams.iIndexDb == 1)  // fragment ion index
   {
      Results* pOutput;

      if (iPrintTargetDecoy != 2)
         pOutput = g_pvQuery.at(iWhichQuery)->_pResults;
      else
         pOutput = g_pvQuery.at(iWhichQuery)->_pDecoys;

      *uiNumTotProteins = (unsigned int)g_pvProteinsList.at(pOutput[iWhichResult].lProteinFilePosition).size();

      int iPrintDuplicateProteinCt = 0; // track # proteins, exit when at iMaxDuplicateProteins

      // get target proteins
      if (iPrintTargetDecoy != 2)  // if not decoy-only
      {
         vector<string> vTmp;      // store decoy matches here to append at end

         comet_fileoffset_t lEntry = pOutput[iWhichResult].lProteinFilePosition;
         for (auto it = g_pvProteinsList.at(lEntry).begin(); it != g_pvProteinsList.at(lEntry).end(); ++it)
         {
            comet_fseek(fpfasta, *it, SEEK_SET);
            if (bReturnFullProteinString)
               fgets(szProteinName, 511, fpfasta);
            else
               fscanf(fpfasta, "%500s", szProteinName);

            szProteinName[500] = '\0';  // limit protein name strings to 500 chars

            if (!strncmp(szProteinName, g_staticParams.szDecoyPrefix, iLenDecoyPrefix))
               vTmp.push_back(szProteinName);
            else
               vProteinTargets.push_back(szProteinName);

            iPrintDuplicateProteinCt++;
            if (iPrintDuplicateProteinCt >= g_staticParams.options.iMaxDuplicateProteins)
               break;
         }

         if (vTmp.size() > 0)      // append any decoy matches now
            vProteinTargets.insert(vProteinTargets.end(), vTmp.begin(), vTmp.end());
      }

      // FIX need to handle decoys
   }
   else  // regular fasta database
   {
      Results *pOutput;

      if (iPrintTargetDecoy != 2)
         pOutput = g_pvQuery.at(iWhichQuery)->_pResults;
      else
         pOutput = g_pvQuery.at(iWhichQuery)->_pDecoys;

      int iPrintDuplicateProteinCt = 0; // track # proteins, exit when at iMaxDuplicateProteins

      *uiNumTotProteins = (unsigned int)(pOutput[iWhichResult].pWhichProtein.size() + pOutput[iWhichResult].pWhichDecoyProtein.size());

      // targets + decoys, targets only, decoys only

      // get target proteins
      if (iPrintTargetDecoy != 2)  // if not decoy-only
      {
         if (pOutput[iWhichResult].pWhichProtein.size() > 0)
         {
            for (auto it=pOutput[iWhichResult].pWhichProtein.begin(); it!=pOutput[iWhichResult].pWhichProtein.end(); ++it)
            {
               comet_fseek(fpfasta, (*it).lWhichProtein, SEEK_SET);
               if (bReturnFullProteinString)
                  fgets(szProteinName, 511, fpfasta);
               else
                  fscanf(fpfasta, "%500s", szProteinName);
               szProteinName[500] = '\0';  // limit protein name strings to 500 chars

               // remove all terminating chars
               while ((szProteinName[strlen(szProteinName) - 1] == '\n') || (szProteinName[strlen(szProteinName) - 1] == '\r'))
                  szProteinName[strlen(szProteinName) - 1] = '\0';

               vProteinTargets.push_back(szProteinName);
               iPrintDuplicateProteinCt++;
               if (iPrintDuplicateProteinCt > g_staticParams.options.iMaxDuplicateProteins)
                  break;
            }
         }
      }

      // get decoy proteins
      if (iPrintTargetDecoy != 1)  // if not target-only
      {
         if (pOutput[iWhichResult].pWhichDecoyProtein.size() > 0)
         {
            // collate decoy proteins, if needed, from target-decoy search
            for (auto it=pOutput[iWhichResult].pWhichDecoyProtein.begin(); it!=pOutput[iWhichResult].pWhichDecoyProtein.end(); ++it)
            {
               if (iPrintDuplicateProteinCt >= g_staticParams.options.iMaxDuplicateProteins)
                  break;

               comet_fseek(fpfasta, (*it).lWhichProtein, SEEK_SET);
               if (bReturnFullProteinString)
                  fgets(szProteinName, 511, fpfasta);
               else
                  fscanf(fpfasta, "%500s", szProteinName);
               szProteinName[500] = '\0';  // limit protein name strings to 500 chars

               // remove all terminating chars
               while ((szProteinName[strlen(szProteinName) - 1] == '\n') || (szProteinName[strlen(szProteinName) - 1] == '\r'))
                  szProteinName[strlen(szProteinName) - 1] = '\0';

               if (strlen(szProteinName) + iLenDecoyPrefix >= WIDTH_REFERENCE)
                  szProteinName[strlen(szProteinName) - iLenDecoyPrefix] = '\0';

               vProteinDecoys.push_back(szProteinName);
               iPrintDuplicateProteinCt++;
               if (iPrintDuplicateProteinCt > g_staticParams.options.iMaxDuplicateProteins)
                  break;
            }
         }
      }
   }
}


// find prev, next AA from first matched protein
// this is only valid if searching indexed db with peptide/protein .idx file
void CometMassSpecUtils::GetPrevNextAA(FILE *fpfasta,
                                       int iWhichQuery,  // which search
                                       int iWhichResult, // which peptide within the search
                                       int iPrintTargetDecoy,    // 0 = target+decoys, 1=target only, 2=decoy only
                                       int iWhichTerm)   // 0=no term constraint, 1=protein N-term, 2=protein C-term
{
   if (g_staticParams.iIndexDb)  // fragment ion or peptide index
   {
      Results *pOutput;

      pOutput = g_pvQuery.at(iWhichQuery)->_pResults;

      pOutput[iWhichResult].cPrevAA = '-';
      pOutput[iWhichResult].cNextAA = '-';

      // for peptide index, if target size is zero, peptide must only be matched to a decoy so don't bother with prev/next AA
      if (g_staticParams.iIndexDb == 2 && g_pvQuery.at(iWhichQuery)->_pResults[iWhichResult].pWhichProtein.size() == 0)
         return;

      if (g_staticParams.iIndexDb == 2)  // peptide index
      {
         for (auto it = pOutput[iWhichResult].pWhichProtein.begin(); it != pOutput[iWhichResult].pWhichProtein.end(); ++it)
         {
            for (auto it2 = g_pvProteinsList.at((*it).lWhichProtein).begin(); it2 != g_pvProteinsList.at((*it).lWhichProtein).end(); ++it2)
            {
               if (SeekPrevNextAA(pOutput, fpfasta, *it2, iWhichQuery, iWhichResult, iWhichTerm))
                  return;
            }
         }

      }
      else // fragment ion index
      {
         comet_fileoffset_t lEntry = pOutput[iWhichResult].lProteinFilePosition;

         for (auto it = g_pvProteinsList.at(lEntry).begin(); it != g_pvProteinsList.at(lEntry).end(); ++it)
         {
            if (SeekPrevNextAA(pOutput, fpfasta, *it, iWhichQuery, iWhichResult, iWhichTerm))
               return;
         }
      }
   }
}


bool CometMassSpecUtils::SeekPrevNextAA(struct Results *pOutput,
                                        FILE *fpfasta,
                                        comet_fileoffset_t tFilePos,
                                        int iWhichQuery,
                                        int iWhichResult,
                                        int iWhichTerm)
{
   string strSeq;
   int iTmpCh = 0;

   int iLenPeptide = (int)strlen(pOutput[iWhichResult].szPeptide);

   comet_fseek(fpfasta, tFilePos, SEEK_SET);

   // skip through protein name string to first carriage return
   while (((iTmpCh = getc(fpfasta)) != '\n') && (iTmpCh != '\r') && (iTmpCh != EOF));

   // Load sequence
   while (((iTmpCh = getc(fpfasta)) != '>') && (iTmpCh != EOF))
   {
      if ('a' <= iTmpCh && iTmpCh <= 'z')
      {
         strSeq += iTmpCh - 32;  // convert toupper case so subtract 32 (i.e. 'A'-'a')
         g_staticParams.databaseInfo.uliTotAACount++;
      }
      else if ('A' <= iTmpCh && iTmpCh <= 'Z')
      {
         strSeq += iTmpCh;
         g_staticParams.databaseInfo.uliTotAACount++;
      }
      else if (iTmpCh == '*')  // stop codon
      {
         strSeq += iTmpCh;
      }

      if (iWhichTerm == 1 && (int)strSeq.length() == iLenPeptide + 1)
      {
         break;  // protein N-terminal peptide identified so have enough sequence now
      }
   }

   if (strSeq.size() < 1)
   {
      printf(" Error: parsed sequence in GetPrevNextAA() is empty.  File pointer %" PRIu64 ", query %d, result %d.\n", tFilePos, iWhichQuery, iWhichResult);
      pOutput[iWhichResult].cPrevAA = pOutput[iWhichResult].cNextAA = '-';
      return false;
   }
   char* szSequence = (char*)malloc(strSeq.size() + 1);

   if (szSequence == NULL)
   {
      printf(" Error: cannot allocate memory for szSequence[%zd]\n", strSeq.size() + 1);
      exit(1);
   }
   strcpy(szSequence, strSeq.c_str());

   int iLenSequence = (int)strlen(szSequence);

   CometSearch cs;
   cs._proteinInfo.iTmpProteinSeqLength = iLenSequence; // used in CheckEnzymeTermini

   if (iWhichTerm == 0)
   {
      size_t iStartPos = 0;

      // Find all occurrences of the peptide in the sequence
      // Take first one consistent with the enzyme constraint for prev/next AA
      while (std::string::npos != (iStartPos = (int)strSeq.find(pOutput[iWhichResult].szPeptide, iStartPos)))
      {
         int iEndPos = (int)iStartPos + iLenPeptide - 1;

         if (cs.CheckEnzymeTermini(szSequence, (int)iStartPos, iEndPos))
         {
            if (iStartPos == 0)
               pOutput[iWhichResult].cPrevAA = '-';
            else
               pOutput[iWhichResult].cPrevAA = szSequence[iStartPos - 1];

            if (iEndPos == iLenSequence - 1)
               pOutput[iWhichResult].cNextAA = '-';
            else
               pOutput[iWhichResult].cNextAA = szSequence[iEndPos + 1];

            free(szSequence);
            strSeq.clear();
            return true;
         }
         else if (g_staticParams.options.bClipNtermMet && iStartPos == 1 && szSequence[0] == 'M' && cs.CheckEnzymeEndTermini(szSequence, iEndPos))
         {
            pOutput[iWhichResult].cPrevAA = 'M';

            if (iEndPos == iLenSequence - 1)
               pOutput[iWhichResult].cNextAA = '-';
            else
               pOutput[iWhichResult].cNextAA = szSequence[iEndPos + 1];

            free(szSequence);
            strSeq.clear();
            return true;
         }

         ++iStartPos;
      }
   }
   else if (iWhichTerm == 1)
   {
      int iEndPos = iLenPeptide; // used for clip n-term met so needs to be set to iLenPeptide

      if (!strncmp(szSequence, pOutput[iWhichResult].szPeptide, iLenPeptide)
         && cs.CheckEnzymeEndTermini(szSequence, iLenPeptide - 1))
      {
         pOutput[iWhichResult].cPrevAA = '-';

         if (iLenSequence >= iLenPeptide)  // for n-term pep, iLenPeptide is following residue position
            pOutput[iWhichResult].cNextAA = szSequence[iLenPeptide];
         else
            pOutput[iWhichResult].cNextAA = '-';

         free(szSequence);
         strSeq.clear();
         return true;
      }
      else if (g_staticParams.options.bClipNtermMet
         && szSequence[0] == 'M'
         && !strncmp(szSequence + 1, pOutput[iWhichResult].szPeptide, iLenPeptide)
         && cs.CheckEnzymeEndTermini(szSequence, iEndPos))
      {
         pOutput[iWhichResult].cPrevAA = 'M';
         pOutput[iWhichResult].cNextAA = szSequence[iEndPos + 1];

         free(szSequence);
         strSeq.clear();
         return true;
      }
   }
   else if (iWhichTerm == 2)
   {
      int iStartPos = iLenSequence - iLenPeptide;

      if (!strncmp(szSequence + iStartPos, pOutput[iWhichResult].szPeptide, iLenPeptide))
      {
         if (cs.CheckEnzymeStartTermini(szSequence, iLenSequence - iLenPeptide))
         {
            if (iStartPos > 0)
               pOutput[iWhichResult].cPrevAA = szSequence[iStartPos - 1];
            else
               pOutput[iWhichResult].cPrevAA = '-';

            pOutput[iWhichResult].cNextAA = '-';

            free(szSequence);
            strSeq.clear();
            return true;
         }
         else if (g_staticParams.options.bClipNtermMet && iStartPos == 1 && szSequence[0] == 'M')
         {
            pOutput[iWhichResult].cPrevAA = 'M';
            pOutput[iWhichResult].cNextAA = '-';

            free(szSequence);
            strSeq.clear();
            return true;
         }
      }
   }

   free(szSequence);
   strSeq.clear();

   return false;
}


// return nth entry in string s
string CometMassSpecUtils::GetField(std::string *s,
                                    unsigned int n,
                                    char cDelimeter)
{
   std::string field;
   std::istringstream isString(*s);

   while ( std::getline(isString, field, cDelimeter) )
   {
      n--;

      if (n == 0)
         break;
   }

   return field;
}


// escape string for XML output
void CometMassSpecUtils::EscapeString(std::string& data)
{
   if (     strchr(data.c_str(), '&')
         || strchr(data.c_str(), '\"')
         || strchr(data.c_str(), '\'')
         || strchr(data.c_str(), '<')
         || strchr(data.c_str(), '>'))
   {
      std::string buffer;
      buffer.reserve(data.size());
      for(size_t pos = 0; pos != data.size(); ++pos)
      {
         switch(data[pos])
         {
            case '&':  buffer.append("&amp;");       break;
            case '\"': buffer.append("&quot;");      break;
            case '\'': buffer.append("&apos;");      break;
            case '<':  buffer.append("&lt;");        break;
            case '>':  buffer.append("&gt;");        break;
            default:   buffer.append(&data[pos], 1); break;
         }
      }
      data.swap(buffer);
   }
}


// input dVal should range from dMin to dMax
char CometMassSpecUtils::NormalizeDoubleToChar(double dVal, double dMin, double dMax)
{
   if (dMax <= dMin)
      return static_cast<char>(-127); // Handle invalid range

   // Normalize dVal to the range [0.0, 1.0] based on dMin and dMax
   double normalizedValue = (dVal - dMin) / (dMax - dMin);

   // Scale to the range [0.0, 255.0]
   double scaledValue = normalizedValue * 255.0;

   // Offset to the char range [-127, 128]
   int charIntValue = static_cast<int>(scaledValue - 127.0);

   // Clamp the value to the valid char range
   if (charIntValue < -127)
      return static_cast<char>(-127);
   else if (charIntValue > 128)
      return static_cast<char>(128);
   else
      return static_cast<char>(charIntValue);
}


double CometMassSpecUtils::DenormalizeCharToDouble(char dChar, double dMin, double dMax)
{
   // Convert the char back to its scaled integer representation [0, 255]
   int scaledIntValue = static_cast<int>(dChar) + 127;

   // Normalize the scaled value to the range [0.0, 1.0]
   double normalizedValue = static_cast<double>(scaledIntValue) / 255.0;

   // Scale back to the original range using dMin and dMax
   double unnormalizedValue = normalizedValue * (dMax - dMin) + dMin;

   return unnormalizedValue;
}