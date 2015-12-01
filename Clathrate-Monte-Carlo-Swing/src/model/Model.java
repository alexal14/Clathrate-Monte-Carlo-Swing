/*****************************************************************************

    Monte Carlo Simulation of sH Clathrate

    Copyright 20014, 2015 Alexander A. Atamas
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/

package model;

public class Model {

	public MetropolisEnsemble simulation = null;
	
	public void setupMonteCarloSimulation(int numMonteCarloRuns, int outputFrequency, double temperature, double pressure, boolean nvt){
		
		Ensemble ensemble;
		
		if ( nvt ){  
			ensemble = Ensemble.NVT;
		} else {
			ensemble = Ensemble.NPT;			
		}
		
		switch(ensemble){
		case NPT:
			simulation = new NPTEnsembleSimulation(numMonteCarloRuns, outputFrequency, temperature, pressure);
			break;
		case NVT:
			simulation = new NVTEnsembleSimulation(numMonteCarloRuns, outputFrequency, temperature);
			break;
		}
	}
	
	public void runMCEngine(int i, int iMol, int iRan){
         
		simulation.engineMonteCarlo(i, iMol, iRan);
	}
	
	
	public Void runMonteCarloSimulation(){
		
		Ensemble ensemble;
		int numMonteCarloRuns;
		
		ensemble = Ensemble.NVT;
		numMonteCarloRuns = 1000000;
		
		switch(ensemble){
		case NPT:
			simulation = new NPTEnsembleSimulation(numMonteCarloRuns);
			break;
		case NVT:
			simulation = new NVTEnsembleSimulation(numMonteCarloRuns);
			break;
		}
		
		if (simulation != null)
			simulation.simRun();
		
		return null;		
	}
	
	public void setGo(boolean flagGo){
	
		simulation.go = flagGo;	
	}
	
	public int getCurrentStepMC(){
		
		return simulation.stepMC;	
	}
	
	public SimulOutputData getSimulOutputData(){
		
		return simulation.simulOutputData;	
	}
	
	public StringBuilder getStringBuilderXYZ(){
		
		return simulation.getStringBuilderAllAtomsXYZ();	
	}

}


