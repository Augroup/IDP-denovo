#include "MurmurHash3.h"
#include "MurmurHash3.cpp"
#include "math.h"
#include <string>
#include <vector>
#include <iostream>

class BloomFilter
{
private:
	int num_of_hashes;
	std::vector<bool> bits;

public:
	BloomFilter(int capacity, double false_positive_probability)
	{
		// SIze of bit array base of the capacity and false_positive_probability
		int bit_array_size = (int)(- (double)(capacity) * log(false_positive_probability)/ ( log(2.0) * log(2.0) ));
		
		// Total num of hashes for the false_positive_probability (research work
		num_of_hashes = (int)((double) bit_array_size * log(2.0) / capacity);
		
		// Init bits all with values false
		bits = std::vector<bool>(bit_array_size, false);
	}

	void add(std::string item)
	{
	    unsigned int result = 0;
	    // Hash each item into and int, and set the corresponding bit
	    // Moding if its over the limit
	    for(int i=0; i<num_of_hashes; i++)
	    {
		// MurmurHash3_x86_32 will return an int representing the hash in result
		MurmurHash3_x86_32(item.c_str(), item.length(), i, &result);
		bits[result%bits.size()] = true;
	    }
	}

	bool contains(std::string item)
	{
	    unsigned int result = 0;
	    for(int i=0; i<num_of_hashes; i++)
	    {
		MurmurHash3_x86_32(item.c_str(), item.length(), i, &result);
		// For an item to be there, all the bits corresponding to the hashes have to be set
		if(!bits[result%bits.size()])
		  return false;
	    }
	    return true;  
	}
};
