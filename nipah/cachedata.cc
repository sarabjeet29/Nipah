#include "cachedata.h"

void cache_Gxy(vector<double>& alpha_vec, vector<double>& nu_vec, double R, int M, int N)
{
    ostringstream dirStringStream, fileStringStream;
    ofstream outputFile;
    const char * filepath;
    const char * dirpath;
    cdouble value;
    
    for(double alpha: alpha_vec)
    {
        for(double nu: nu_vec)
        {
            cout << "alpha: " << alpha <<", nu: "<<nu<<"\n";
            dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
            dirpath = dirStringStream.str().c_str();
            mkdir(dirpath,0777);
            
            fileStringStream << dirStringStream.str().c_str() << "/Gxy_M=" << M << "_N=" << N << ".txt";
            filepath = fileStringStream.str().c_str();
            outputFile.open( filepath );
            if( abs(round(nu / alpha ) - nu / alpha) < 1e-12 )
            {
                cout <<"adjusted for alpha: "<<alpha<<", nu: "<<nu<< "\n";
                nu += 1e-6;
            };
            for(int i = 0; i < M; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    value = Gxy(cdouble(R,0) * exp(cdouble(0,2*PI*i/M)),
                                cdouble(R,0) * exp(cdouble(0,2*PI*j/N)), alpha, nu);
                    outputFile << real(value) << ' ' << imag(value) <<'\n';
                };
            };
            outputFile.close();
            dirStringStream.str(string());
            fileStringStream.str(string());
        };
    };
};

void save_prob_alpha_nu(vector<double>& alpha_vec, vector<double>& nu_vec,
                        vector<int>& outbsPrimaryInfo, vector<int>& outbsNoPrimaryInfo,
                        int sumPrimaryOutbs, int M, int N)
{
    double pij;
    vector<cdouble> cachedGxy;
    ofstream outputFile;
    outputFile.open("./results/temp_prob_alpha_nu.txt");
    ostringstream pijstr;
    for(double alpha : alpha_vec)
    {
        vector<double> probRow;
        for( double nu : nu_vec)
        {
            cout << "alpha: " << alpha <<", nu: "<< nu <<"\n";
            read_cached_Gxy(alpha, nu, M, N, cachedGxy);

            pij = getProbObsOutbs(alpha, nu, outbsPrimaryInfo, outbsNoPrimaryInfo, sumPrimaryOutbs, N, M, cachedGxy);
            cachedGxy.clear();
            outputFile << pij << ' ';
        };
        outputFile << '\n';
    };
    outputFile.close();
};

void save_PmnPrimaryInfo(vector<double> alpha_vec, vector<double> nu_vec,
                         vector<int> outbsPrimaryInfo, int sumPrimaryOutbs,
                         int M, int N)
{
    ostringstream dirStringStream, fileStringStream;
    ofstream myfile;
    vector<double> Pmn;
    const char * filepath;
    const char * dirpath;
    vector<cdouble> cachedGxy;
    
    for(double alpha: alpha_vec)
    {
        for(double nu: nu_vec)
        {

            dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
            dirpath = dirStringStream.str().c_str();
            mkdir(dirpath,0777);
            
            dirStringStream << "/p";
            dirpath = dirStringStream.str().c_str();
            mkdir(dirpath,0777);
            
            cout << "alpha: " << alpha <<", nu: "<<nu<<"\n";
            read_cached_Gxy(alpha, nu, M, N, cachedGxy);
            for(int s: outbsPrimaryInfo)
            {
                fileStringStream << dirStringStream.str().c_str() << "/" << s <<".txt";
                filepath = fileStringStream.str().c_str();
                fileStringStream.str(string());
                myfile.open ( filepath );
                getPmPrimaryInfo(outbsPrimaryInfo, sumPrimaryOutbs, s, N, M, cachedGxy, Pmn);

                for( double pij : Pmn)
                    myfile << pij <<'\n';
                Pmn.clear();
                myfile.close();
            };
            Pmn.clear();
            cachedGxy.clear();
            dirStringStream.str(string());
        };
    };
};

void save_PmnNoPrimaryInfo(vector<double> alpha_vec, vector<double> nu_vec,
                           vector<int> outbsNoPrimaryInfo, int M, int N)
{
    ostringstream dirStringStream, fileStringStream;
    ofstream myfile;
    const char * filepath;
    const char * dirpath;
    vector<vector<double>> Pmn;
    vector<cdouble> cachedGxy;
    
    for(double alpha: alpha_vec)
    {
        for(double nu: nu_vec)
        {
            dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
            dirpath = dirStringStream.str().c_str();
            mkdir(dirpath,0777);
            
            dirStringStream << "/np";
            dirpath = dirStringStream.str().c_str();
            mkdir(dirpath,0777);
            
            cout << "alpha: " << alpha <<", nu: "<<nu<<"\n";
            read_cached_Gxy(alpha, nu, N, N, cachedGxy);
            getPmNoPrimaryInfo(M, N, cachedGxy, Pmn);
            cachedGxy.clear();
            for( int s: outbsNoPrimaryInfo)
            {
                fileStringStream << dirStringStream.str().c_str() << "/" << s <<".txt";
                filepath = fileStringStream.str().c_str();
                fileStringStream.str(string());
                myfile.open ( filepath );
                
                for( double pij : Pmn[s])
                    myfile << pij <<'\n';
                
                myfile.close();
            };
            Pmn.clear();
            dirStringStream.str(string());
        };
    };
};

void cache_Hxzuw(vector<double>& alpha_vec, vector<double>& nu_vec, vector<double>& p_vec,
                 vector<double> q_vec, double R, double M, double N)
{
    ostringstream dirStringStream, fileStringStream;
    gzFile outputFile;
    const char * filepath;
    const char * dirpath;
    cdouble value;
    ostringstream strvalue;
    strvalue << " ";
    for(double alpha: alpha_vec)
    {
        for(double nu: nu_vec)
        {
            cout << "alpha: " << alpha << ", nu: "<< nu << "\n" ;
            dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
            dirpath = dirStringStream.str().c_str();
            mkdir(dirpath,0777);

            for(double p: p_vec)
            {
                for(double q: q_vec)
                {
                    cout <<"\tp: " <<p<< ", q: " <<q<<"\n";
                    dirStringStream << "/p=" << int(p*100)/100.0 << ",q=" << int(q*100)/100.0;
                    dirpath = dirStringStream.str().c_str();
                    mkdir(dirpath,0777);
                    
                    fileStringStream << dirStringStream.str().c_str() << "/Hxzuw_M=" << M <<
                    "_N=" << N << ".gz";
                    filepath = fileStringStream.str().c_str();
                    outputFile = gzopen(filepath, "wb9");
                    if( abs(round(nu / alpha ) - nu / alpha) < 1e-12 )
                    {
                        cout <<"adjusted for alpha: "<<alpha<<", nu: "<<nu<< "\n";
                        nu += 1e-6;
                    };
                    for(int i = 0; i < M; i++)
                    {
                        for(int j = 0; j < M; j++)
                        {
                            for(int k = 0; k < N; k++)
                            {
                                for(int l = 0; l < N; l++)
                                {
                                    value = Hxzuw(cdouble(R,0) * exp(cdouble(0,2*PI*j/M)),
                                                  cdouble(R,0) * exp(cdouble(0,2*PI*l/M)),
                                                  cdouble(R,0) * exp(cdouble(0,2*PI*i/N)),
                                                  cdouble(R,0) * exp(cdouble(0,2*PI*k/N)),
                                                  alpha, nu, p, q);
                                    strvalue << real(value) << ' ' << imag(value) << "\n";
                                    gzprintf(outputFile, strvalue.str().c_str());
                                    strvalue.str(string());
                                };
                            };
                        };
                    };
                    gzclose(outputFile);
                    dirStringStream.str(string());
                    fileStringStream.str(string());
                };
            };
        };
    };
};

//void save_prob_p_q(vector<double>& alpha_vec, vector<double>& nu_vec,
//                   vector<double>& p_vec, vector<double>& q_vec,
//                   vector<int>& outbsPrimaryInfo, vector<int>& outbsNoPrimaryInfo,
//                   vector<int>& deathsPrimaryInfo, vector<int>& deathsNoPrimaryInfo,
//                   int sumPrimaryOutbs, int sumPrimaryDeaths, int M, int N, double R )
//{
//    ostringstream dirStringStream, fileStringStream;
//    ofstream myfile;
//    const char * filepath;
//    const char * dirpath;
//    vector<vector<double>> probObsDeaths;
//    vector<double> probRow;
//    double pij;
//    for(double alpha: alpha_vec)
//    {
//        for( double nu: nu_vec)
//        {
//            dirStringStream << "./results/alpha=" << int(alpha*100)/100.0 << ",nu=" << int(nu*100)/100.0;
//            dirpath = dirStringStream.str().c_str();
//            mkdir(dirpath,0777);
//            
//            fileStringStream << dirStringStream.str().c_str() <<"/prob_p_q.txt";
//            filepath = fileStringStream.str().c_str();
//            fileStringStream.str(string());
//            dirStringStream.str(string());
//            myfile.open ( filepath );
//            
//            if( abs(round(nu / alpha ) - nu / alpha) < 1e-12 )
//            {
//                cout << "slightly adjusted nu for alpha: "<<alpha<<", nu: "<<nu<< "\n";
//                nu += 1e-6;
//            };
//            
//            for(double p: p_vec)
//            {
//                for(double q: q_vec)
//                {
//                    cout << "alpha: " << alpha <<", nu: "<< nu << ", p: " << p <<", q: "<< q << "\n";
//                    
//                    pij = getProbObsOutbsAndDeaths(outbsPrimaryInfo, outbsNoPrimaryInfo,
//                                                   deathsPrimaryInfo, deathsNoPrimaryInfo,sumPrimaryOutbs, sumPrimaryDeaths, M, N,
//                                                   cachedHxzuw);
//                    myfile << pij << ' ';
//                    probRow.push_back(pij);
//                };
//                probObsDeaths.push_back(probRow);
//                probRow.clear();
//                myfile << '\n';
//            };
//            probObsDeaths.clear();
//        };
//    };
//    myfile.close();
//};
