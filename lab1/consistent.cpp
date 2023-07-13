#include <iostream>
#include <vector>

using namespace std;
typedef long long ll;

ll consistent(int &argc, char **&argv) {
    vector<int> a, b;
    ll n = atoi(argv[1]);
    for (ll i = 0; i < n; ++i) {
        a.push_back(1);
        b.push_back(1);
    }
    ll s=0;
    for(ll i=0;i<n;i++)
        for(ll l=0;l<n;l++)
            s+=a[i]*b[l];
    return s;
}

int main(int argc, char **argv) {
    cout<<consistent(argc,argv)<<endl;
}