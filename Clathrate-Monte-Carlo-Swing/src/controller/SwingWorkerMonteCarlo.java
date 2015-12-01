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

package controller;

import java.awt.geom.Point2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.SwingWorker;

import model.AtomVector;
import model.IMolecularSimulation;
import model.Model;
import model.SimulOutputData;
import view.View;

public class SwingWorkerMonteCarlo extends SwingWorker<Void, SimulOutputData> implements PropertyChangeListener, IMolecularSimulation {

	private Model theModel;
	private View  theView;
	private Controller  theController;
	private int numbMCSimulations;
	private int outputFrequency;
	private double temperature;
	private double pressure;
	private boolean nvt;
	private StringBuilder sbAllAtomsXYZ = new StringBuilder();
	
	private int numbOfMolsActual;
	
	SwingWorkerMonteCarlo(Controller theController, Model theModel, View theView, int numbMCSimulations, int outputFrequency, double temperature, double pressure, boolean nvt){
		
		this.theModel = theModel;
		this.theView  = theView;
		this.theController = theController;
		this.numbMCSimulations = numbMCSimulations;
		this.outputFrequency = outputFrequency;
		this.temperature = temperature;
		this.pressure = pressure;
		this.nvt = nvt;
		
		addPropertyChangeListener(this);
	}
	
	@Override
	protected Void doInBackground() throws Exception {
		
		//theModel.runMonteCarloSimulation();
		
		numbOfMolsActual =  272;
		
		theModel.setupMonteCarloSimulation(numbMCSimulations, outputFrequency, temperature, pressure, nvt);
		
		int iMol, iRan;
		int start, end;
		Random random = new Random();
		
		setProgress(0);
		
		for (int i=0; i < numbMCSimulations; i++){
			
			if ( isCancelled() ) break;
		
			start = 0;
		    end = numbOfMolsActual - 1;		    
		    iMol = RandomInteger(start, end, random);   // generate random index of molecule
		    
			start = 1;
		    end = 2*numbOfMolsActual + 3;		    
		    iRan = RandomInteger(start, end, random);
		    
		    theModel.runMCEngine(i, iMol, iRan);
		    
		    int istep = theModel.getSimulOutputData().i;
		    double en = theModel.getSimulOutputData().energy;
		    double enH = theModel.getSimulOutputData().energyH; 
		    AtomVector box = theModel.getSimulOutputData().box;
		    
			if ((i % outputFrequency) == 0) {
				System.out.printf("i,E,box = %9d    %8.2f   %8.2f   %s\n", istep, en , enH, box);
			}
			
			if ( i % outputFrequency == 0 ){
				
//				System.out.println("publish");
				setProgress(100 * i/numbMCSimulations);
				publish(theModel.getSimulOutputData());
			}   
		}	
		
		return null;
	}
	
	@Override	
	protected void process(List<SimulOutputData> data){
		
		SimulOutputData simulOutput = data.get(data.size() - 1);
		
		int i = simulOutput.i;
		if (i==0) i++;
		
	    double enU = theModel.getSimulOutputData().energy;
	    double enH = theModel.getSimulOutputData().energyH;
	    double xBox = theModel.getSimulOutputData().box.x;
	    double yBox = theModel.getSimulOutputData().box.y;
	    double zBox = theModel.getSimulOutputData().box.z;
	    ArrayList<Point2D> gr = theModel.getSimulOutputData().gr;
	    sbAllAtomsXYZ = theModel.getSimulOutputData().sbAllAtomsXYZ;
	    
	    theController.gr = gr;
	    theController.sbAllAtomsXYZ = sbAllAtomsXYZ;
	    
	    double vol = xBox*yBox*zBox;
	    theController.nOverVolume = (double)numbOfMolsActual/vol;
	    double Mw = 18.01528;    					// Mw, molar mass of water molecule, g/mol 
	    Mw = Mw*1e-3;            					// molar mass of water in kg/mol
	    double Mass = numbOfMolsActual*Mw;
	    double rho = Mass/vol*1e30/NA/1e3; 			// density, g/cm^3	
	    
		theView.setSimulProgressText(i, numbMCSimulations);
	    theView.updateGUI(i, enU, enH, xBox, yBox, zBox, gr, rho, sbAllAtomsXYZ);
		
	}
	
	@Override
	protected void done(){
		
		if ( isCancelled() ){
			
			theView.setSimulCancelledText();
			
		} else {
			
			theView.setSimulCompletedText();	

		}	
		
		theView.getProgressBarSimulation().setValue(theView.getProgressBarSimulation().getMaximum());
	}
	
	
	private static int RandomInteger(int aStart, int aEnd, Random aRandom) {
		if (aStart > aEnd) {
			throw new IllegalArgumentException("Start cannot exceed End.");
		}
		// get the range, casting to long to avoid overflow problems
		long range = (long) aEnd - (long) aStart + 1;
		// compute a fraction of the range, 0 <= frac < range
		long fraction = (long) (range * aRandom.nextDouble());
		int randomNumber = (int) (fraction + aStart);
		return randomNumber;
	}

	@Override
	public void propertyChange(PropertyChangeEvent evt) {

		if ("progress" == evt.getPropertyName()){
		
			int progress = (Integer) evt.getNewValue();
			theView.getProgressBarSimulation().setValue(progress);
		}		
		
	}

}
