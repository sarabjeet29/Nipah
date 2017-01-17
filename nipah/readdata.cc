#include "readdata.h"

void read_cached_Gxy(double alpha, double nu, int M, int N, vector<cdouble>& Gxy)
{
    ostringstream dirStringStream, fileStringStream;
    ifstream inputFile;
    const char * filepath;
    const char * dirpath;
    string line;
    double realval, imagval;

    dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
    dirpath = dirStringStream.str().c_str();
    fileStringStream << dirStringStream.str().c_str() << "/Gxy_M=" << M << "_N=" << N << ".txt";
    filepath = fileStringStream.str().c_str();
    inputFile.open( filepath );
    int i = 0;
    if (inputFile.is_open())
    {
        while ( getline (inputFile, line) )
        {
            istringstream ss(line);
            ss >> realval >> imagval;
            Gxy.push_back(cdouble(realval, imagval));
            i++;
        }
        inputFile.close();
    }
}