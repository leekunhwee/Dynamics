#include<stdlib.h>
#include<iostream>
#include<fstream>

using namespace std;

/***************************************************************************/

int GetWordsInLine(istream &stream, int &nbrmotslus, char **Mot)

/* Fonction permettant d'extraire d'une ligne tous les mots separes par un
   ou plusieurs espaces et tabulations */

  {
  char c;
  int lncrtword=0; // longueur du mot courant;
  nbrmotslus=0; // nombre de mots
  if (stream.rdstate()) return 1;
  do
    {
    // Reading spaces
    while ((stream.peek()==' ') || (stream.peek()=='\t')) c=stream.get();
    // Reading word
    nbrmotslus++; lncrtword=0;
    while ((stream.peek()!=' ') && (stream.peek()!='\t') &&
           (stream.peek()!='\n') && (stream.peek()!='\r') &&
           (stream.peek()!='\f') && (stream.peek()!=EOF))
      {
      Mot[nbrmotslus-1][lncrtword]=stream.get();
      lncrtword++;
      }
    if (lncrtword==0) nbrmotslus--;
    else Mot[nbrmotslus-1][lncrtword]='\0';
    }
  while ((stream.peek()!='\n') && (stream.peek()!='\f') 
         &&(stream.peek()!='\r') && !stream.rdstate()); 
  // Reading end of line characters
  while ((stream.peek()=='\n') || (stream.peek()=='\f') 
         || (stream.peek()=='\r')) c=stream.get(); 
  // Looking for end of file before returning
  if (stream.peek()==EOF) c=stream.get(); 
  return 0;
  } // fin de GetWords

/***************************************************************************/

int FindNextMinusOne(istream &stream, char **Mot)
{
int NbrWords;
do
  {
  if (GetWordsInLine(stream,NbrWords,Mot)) return 1;
  }
  while ((NbrWords!=1) || (strcmp(Mot[0],"-1")!=0));
return 0;
}

/***************************************************************************/

int FindNextUff(istream &stream, int type, char **Mot)
{
int NbrWords;
 do
   {
   if (FindNextMinusOne(stream,Mot)) return 1;
   if (GetWordsInLine(stream,NbrWords,Mot)) return 1;
   if ((NbrWords==1) && (atoi(Mot[0])==type))
     return 0;
   else if (FindNextMinusOne(stream,Mot)) return 1;
   }
   while (!stream.rdstate());
 return 1;
}

/***************************************************************************/

main(int argc, char *argv[])

{

if (argc!=2)
  {
  cout << "Usage: decode file" << endl;
  cout << "Affiche les caracteres hexadecimaux" << endl;
  exit(1);
  }

ifstream infile(argv[1]);
if (!infile)
  {
  cout << "Cannot open file " << argv[1] << endl;
  exit(1);
  }


int i;
 unsigned char c; unsigned int code; int nbrmots;
 char **Mot;
 Mot =new char*[20];
 for (i=0; i<20; i++) Mot[i]=new char[30];

while (!infile.eof())
  {
    GetWordsInLine(infile,nbrmots,Mot);
    for (i=0;i<nbrmots;i++) cout << " " << Mot[i];
    cout << endl;
  }

infile.close();

return 0;

}
