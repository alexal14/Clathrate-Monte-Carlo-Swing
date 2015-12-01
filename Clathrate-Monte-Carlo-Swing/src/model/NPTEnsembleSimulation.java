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

public class NPTEnsembleSimulation extends MetropolisEnsemble {

	private int outputFrequency;
	
	NPTEnsembleSimulation(final int numbMCSimulations){
		super(numbMCSimulations);
		
		this.outputFrequency = 2500;
	}
	
	NPTEnsembleSimulation(final int numbMCSimulations, final int outputFrequency, final double temperature, final double pressure){
		super(numbMCSimulations, temperature, pressure);
		
		this.outputFrequency = outputFrequency;
	}
	
	@Override
	protected void engineMonteCarlo(int i, int iMol, int iRan) {
		
	    if(iRan > 2*numbOfMolsActual){
	    	
		       enTot = volumeChange(molModelMap, 
		    		  			    molecules, 
	                                numbOfMolsActual,
	                                box,
	                                rcutSq,
	                                tableEnergyAllMoleciles,
	                                enTot
		    		  			   );
		    	
		    } else if (iRan < numbOfMolsActual){
		    	
				enTot = moleculeTranslation(iMol,
				 							molModelMap, 
				 							molecules, 
				 							numbOfMolsActual,
				 							box,
				 							rcutSq,
				 							tableEnergyAllMoleciles,
				 							enTot
										   );
		    	
		    } else {
		    	
				enTot = moleculeRotation(iMol,
		    			 				 molModelMap, 
		    			 				 molecules, 
		    			 				 numbOfMolsActual,
		    			 				 box,
		    			 				 rcutSq,
		    			 				 tableEnergyAllMoleciles,
		    			 				 enTot
		   								);		    	
		    }
	    
	    rdf();
	    
    
/*	    if ((stepMC % freqOfSysOut) == 0) {
	    	System.out.printf("i,E,box = %9d    %8.2f    %8.3f    %8.3f    %8.3f \n", stepMC, enTot / (double) numbOfMolsActual, box.x, box.y, box.z);
	    	writeRDF_File("output_simulation_g_Ow_Hw.xyz", box, numbWaterMols, stepMC, g);
		}*/
	    
		if ((i % outputFrequency) == 0) {
			double vol = box.x * box.y * box.z;
		    double press = 1.0;
		    double pv = (1e-3)*press*(1e8)*NA*vol*(1e-30)*0.2388458966275 ; 
			double enH = (enTot + pv) / (double) numbOfMolsActual;
			
			writeRDF_File("output_simulation_g_Ow_Hw_NPT.xyz", box, numbWaterMols, i, g);
			normalizeRDF(box, numbWaterMols, i, g);		
			
			StringBuilder sb = getStringBuilderAllAtomsXYZ();
			
			simulOutputData.setData(i, enTot / (double) numbOfMolsActual, enH, box, gr, sb);		
	
		}
		
	}
	
	protected void simulationRun(){
		
		simRun();		
	}
	

}
