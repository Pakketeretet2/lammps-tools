#ifndef DUMP_INTERPRETER_DCD_H
#define DUMP_INTERPRETER_DCD_H

#include "block_data.h"
#include "dump_reader.h"
#include "reader_core_binary.h"

class dump_interpreter_dcd : public dump_interpreter
{
public:
	// To convert NAMD units to A/ps/fs.
	static constexpr double PDBVELVACTOR = 20.45482706;
	static constexpr double TIMEFACTOR   = 48.88821;
	

	dump_interpreter_dcd() : namd_units(true)
	{}
	virtual ~dump_interpreter_dcd(){}

	
	virtual int next_block( reader_core *r, block_data &block );
	
        
	
	
	
private:
	bool namd_units;
};






#endif // DUMP_INTERPRETER_DCD_H
