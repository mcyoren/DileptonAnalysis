////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// See WriteEvent.h for detailed discription 
//
// used to create particle list so as to write ROOT and Oscar output 
//
// Roli Esha 05/13/2022
//

#include "WriteEvent.h"
#include <fstream>
#include <iostream>


void WriteEvent::ClearEvent(){

	stable_particles = -999;

	WriteTrackList.clear();

}

std::ofstream* OpenOutputFile(const char * file){

  	std::cout << "Opening output file: " << file << std::endl;
  	std::cout << std::endl;

  	std::ofstream * output_file = new std::ofstream;;
  	output_file->open(file);

  	return (output_file);
}

void WriteOscarEvent(std::ofstream* file, WriteEvent Event){

	int parentID = -999;
	WriteTrack Track;

	*file << 0 << "\t" << Event.GetNEntries() << std::endl;

	for(int i = 0; i < Event.GetNEntries(); i++){

		Track = Event.GetWriteTrack(i);

		if(Track.GetFinal() == 0){

			parentID = i;

			*file << Track.GetNum() << "\t" << "0" << "\t" << Track.GetID() << "\t" << Track.GetPx() << "\t" << Track.GetPy() << "\t" << Track.GetPz() << "\t" << Track.GetEnergy() << "\t" << Track.GetMass() << "\t" << Track.GetXpos() << "\t" << Track.GetYpos() << "\t" << Track.GetZpos() << "\t" << Track.GetBranch() << std::endl;
		}

		else if(Track.GetFinal() == 1){

			*file << Track.GetNum() << "\t" << Track.GetID() << "\t" << parentID << "\t" << Track.GetPx() << "\t" << Track.GetPy() << "\t" << Track.GetPz() << "\t" << Track.GetEnergy() << "\t" << Track.GetMass() << "\t" << Track.GetXpos() << "\t" << Track.GetYpos() << "\t" << Track.GetZpos() << "\t" << "0" << std::endl;

		}
	}

	*file << "0" << "\t" << "0" << std::endl;

}



