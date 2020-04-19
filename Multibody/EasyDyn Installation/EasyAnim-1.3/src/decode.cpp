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

  while ((stream.peek()=='\n') || (stream.peek()=='\f') 
         && (stream.peek()=='\r')) c=stream.get(); 
  if (stream.peek()==EOF) c=stream.get();
  } // fin de GetWords

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
unsigned char c; unsigned int code;
while (!infile.eof())
  {
  c=infile.get(); code=c;
  if (!infile.eof())
    {
    cout << " " << code;
    }
  }

infile.close();

return 0;

}
