#include <iostream>
#include <fstream>
using namespace std;

void XoutofNg()
{
	int x = 6;
	int n = 8;


	if( x == 6 && n == 10 )
	{
	// calculates all possible configurations of 6 from 10 available options. This is to be used to 
	// separate the eta--> 3pi0 --> 6g from 10 detected g (eta' --> eta pi0 pi0 where the latter two pi0 decays to 4g) 
		ofstream comb;
		comb.open("combsixoften.txt");

		for( int a = 0; a < n; a++ )
		{
			for( int b = ( a + 1 ); b < n; b++ )
			{
				for( int c = ( b + 1 ); c < n; c++ )
				{
					for( int d = ( c + 1 ); d < n; d++ )
					{
						for( int e = ( d + 1 ); e < n; e++ )
						{
							for( int f = ( e + 1 ); f < n; f++ )
							{
								if( f < 9 )
									comb << "{ " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << " }" << ", \t";
								else
									comb << "{ " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << " }" << ", \n";			
							}
						}
					}
				}
			}
		}
	}


	if( x == 6 && n == 7 )
	{
	// calculates all possible configurations of 6 from 10 available options. This is to be used to 
	// separate the eta--> 3pi0 --> 6g from 10 detected g (eta' --> eta pi0 pi0 where the latter two pi0 decays to 4g) 
		ofstream comb;
		comb.open("combsixofseven.txt");

		for( int a = 0; a < n; a++ )
		{
			for( int b = ( a + 1 ); b < n; b++ )
			{
				for( int c = ( b + 1 ); c < n; c++ )
				{
					for( int d = ( c + 1 ); d < n; d++ )
					{
						for( int e = ( d + 1 ); e < n; e++ )
						{
							for( int f = ( e + 1 ); f < n; f++ )
							{
									comb << "{ " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << " }" << ", \n";			
							}
						}
					}
				}
			}
		}
	}
	
	
		if( x == 6 && n == 8 )
	{
	// calculates all possible configurations of 6 from 10 available options. This is to be used to 
	// separate the eta--> 3pi0 --> 6g from 10 detected g (eta' --> eta pi0 pi0 where the latter two pi0 decays to 4g) 
		ofstream comb;
		comb.open("combsixofeight.txt");

		for( int a = 0; a < n; a++ )
		{
			for( int b = ( a + 1 ); b < n; b++ )
			{
				for( int c = ( b + 1 ); c < n; c++ )
				{
					for( int d = ( c + 1 ); d < n; d++ )
					{
						for( int e = ( d + 1 ); e < n; e++ )
						{
							for( int f = ( e + 1 ); f < n; f++ )
							{
									comb << "{ " << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << " }" << ", \n";			
							}
						}
					}
				}
			}
		}
	}

	comb.close();
}

