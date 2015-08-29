//
// Created by feixh on 8/29/15.
//

#include <fstream>
#include <iostream>
#include <vector>

#include <pnpsolvers/P3p.h>

using namespace std;

int main()
{
    ifstream ifs("../test/data.txt", ifstream::in );
    assert( ifs.is_open() );
    vector< float > pc;
    vector< float > pts;
    int n, n_inliers;
    ifs >> n >> n_inliers;
    for ( int i = 0; i < n; ++i )
    {
        float x, y;
        ifs >> x >> y;
        if ( i < n_inliers ){
            pc.push_back( x );
            pc.push_back( y );
        }
    }
    for ( int i = 0; i < n; ++i )
    {
        float x, y, z;
        ifs >> x >> y >> z;
        if ( i < n_inliers ){
            pts.push_back( x );
            pts.push_back( y );
            pts.push_back( z );
        }
    }
    float gwc[3][4];
    for ( int i = 0; i < 3; ++i ) {
        for (int j = 0; j < 4; ++j) {
            ifs >> gwc[i][j];
        }
    }
    ifs.close();
}