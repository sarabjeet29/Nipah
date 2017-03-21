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

void read_cached_Hxzuw(double alpha, double nu, double p, double q, int M, int N,
                       vector<cdouble>& Hxzuw)
{
    ostringstream dirStringStream, fileStringStream;
    gzFile inputFile;
    const char * filepath;
    const char * dirpath;
    double realval, imagval;
    string line;

    dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0 << "/p=" << int(p*100)/100.0 << ",q=" << int(q*100)/100.0;;
    dirpath = dirStringStream.str().c_str();
    fileStringStream << dirStringStream.str().c_str() << "/Hxzuw_M=" << M << "_N=" << N << ".gz";
    filepath = fileStringStream.str().c_str();
    inputFile = gzopen(filepath, "rb");
    
    unsigned int size;
    gzread(inputFile, (void*) &size, sizeof(size));
    string data;
    data.resize(size / sizeof(char));
    gzread(inputFile, (void*) data.data(), size);
    
    istringstream iss(data);
    string token;
    while(getline(iss, token, '\n'))
    {
        istringstream values(token);
        values >> realval >> imagval;
        Hxzuw.push_back(cdouble(realval, imagval));
    }
    Hxzuw[0] = cdouble(1,0);
    gzclose(inputFile);
}