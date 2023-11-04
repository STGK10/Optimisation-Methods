/*
4 4  # Number of suppliers and consumers
30 40 25 35  # Supply values for suppliers
20 30 25 35  # Demand values for consumers
4.0 5.0 6.0 8.0  # Cost matrix (cost per unit transportation)
7.0 6.0 9.0 7.0
3.0 8.0 5.0 6.0
9.0 7.0 4.0 5.0


 */
#include <iostream>
#include "TransportTask.h"


int main(void)
{
	TransportTask tt;
	
	if (tt.generate("task1.txt"))
	{
	    tt.run();
	}
	else
	{
		std::cerr << "Something gona wrong!" << std::endl;
	}
	
	return 0;
}
