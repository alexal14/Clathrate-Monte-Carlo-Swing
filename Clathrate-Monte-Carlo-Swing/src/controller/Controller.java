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

import model.*;
import view.*;

import java.awt.CardLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import javax.swing.JOptionPane;
import javax.swing.JToggleButton;

public class Controller implements ActionListener, ItemListener {

	private Model theModel;
	private View  theView;
	
	public ArrayList<Point2D> gr = new ArrayList<Point2D>();
	public ArrayList<Point2D> grDensity = new ArrayList<Point2D>();
	public double nOverVolume;
	public StringBuilder sbAllAtomsXYZ = new StringBuilder();
	private SwingWorkerMonteCarlo swingWorker = null;
	
	public Controller(Model theModel, View theView) {
		
		this.theModel = theModel;
		this.theView  = theView;
		
		theView.getBtnRun().addActionListener(this);
		theView.getBtnStop().addActionListener(this);
		theView.getMntmQuit().addActionListener(this);
		theView.getMntmAbout().addActionListener(this);
		theView.getRdbtnNVT().addActionListener(this);
		theView.getRdbtnNPT().addActionListener(this);
		theView.getBtnRDF().addActionListener(this);
		theView.getBtnNormedRDF().addActionListener(this);
		theView.getTglBtnBox().addItemListener	(this);
		theView.getTglBtnPerspective().addItemListener(this);
		
		theView.showTempPressFreqNPT();
		
////////////////////////////////////////////////////////////////////////////////////////////////////////

		theView.getBtnRDF().setEnabled(false);
		theView.getBtnNormedRDF().setEnabled(false);
		theView.emptyOutputTable();

		int numbMC = theView.getNumbMCSteps();
		int outputFrequency = theView.getFrequency();
		double temperature = theView.getTemperature();
		double pressure = theView.getPressure();

		boolean nvt = false;
		if (theView.isNVTSelected()) {
			nvt = true;
		} else if (theView.isNPTSelected()) {
			nvt = false;
		}

		theModel.setupMonteCarloSimulation(numbMC, outputFrequency, temperature, pressure, nvt);
		
		sbAllAtomsXYZ.setLength(0);
		sbAllAtomsXYZ = theModel.getStringBuilderXYZ();		
		theView.showCrystalInJMol(sbAllAtomsXYZ);

//////////////////////////////////////////////////////////////////////////////////////////////////////
		
	}
	
	public void actionPerformed(ActionEvent ae) {

		Object btnClicked = ae.getSource();
		
		if ( btnClicked.equals(theView.getBtnRun()) ){
			
			theView.getBtnRun().setEnabled(false);
			theView.getBtnStop().setEnabled(true);	
			theView.getBtnRDF().setEnabled(true);
			theView.getBtnNormedRDF().setEnabled(true);
			theView.emptyOutputTable();
			
			int numbMC = theView.getNumbMCSteps();			
			int outputFrequency = theView.getFrequency();
			double temperature = theView.getTemperature();
			double pressure = theView.getPressure();
			
			boolean nvt = false;			
			if ( theView.isNVTSelected() ){
				nvt = true;
			} else if ( theView.isNPTSelected() ){
				nvt = false;
			}
			
			swingWorker = new SwingWorkerMonteCarlo(this, theModel, theView, numbMC, outputFrequency, temperature, pressure, nvt);
			swingWorker.execute();
			
			return;
			
		} else if ( btnClicked.equals(theView.getBtnStop()) ) {
			
			theView.getBtnRun().setEnabled(true);
			theView.getBtnStop().setEnabled(false);
			
			//theModel.setGo(false);
			swingWorker.cancel(true);
			
			return;
			
		} else if( btnClicked.equals(theView.getMntmQuit()) ) {
			
			theView.quitApp();
			return;
			
        } else if( btnClicked.equals(theView.getMntmAbout()) ) {
        	
//        	JOptionPane.showMessageDialog(theView, "Button \"About\" was clicked");
        	theView.showAboutDialog();
			return;
			
        } else if( btnClicked.equals(theView.getRdbtnNVT()) ) {        	

        	theView.showTempFreqNVT();

			return;
			
        } else if( btnClicked.equals(theView.getRdbtnNPT()) ) {
        	
        	theView.showTempPressFreqNPT();
        	
			return;
			
        } else if( btnClicked.equals(theView.getBtnRDF()) ) {
        	
        	if ( swingWorker != null ){
        		theView.showGraphsDialog(gr, "Radial Distribution Function (RDF)", "g(r) O-H");
        	}
        	
			return;
			
        } else if( btnClicked.equals(theView.getBtnNormedRDF()) ) {
        	
        	if ( swingWorker != null ){
        		
        		grDensity.clear();
        		
        		double grNorm, r, r2;
        		
                for (int j = 0; j < gr.size(); j++) {
                	
                	r = gr.get(j).getX();
                	r2 = r*r;
                	grNorm = gr.get(j).getY()*4*Math.PI*r2*nOverVolume;
                	
                	Point2D pt = new Point2D.Double(r, grNorm);   
                	
                	grDensity.add(pt);
        		} 
        		
        		theView.showGraphsDialog(grDensity, "Particle distribution", "g(r) O-H normalized ");
        	}
        	
			return;
			
        } 		
	}

	@Override
	public void itemStateChanged(ItemEvent e) {

		JToggleButton tBtn = (JToggleButton)e.getSource();
		
		if( tBtn.equals(theView.getTglBtnBox()) ) {      	

		      if(e.getStateChange()==ItemEvent.SELECTED){
		          
		          theView.setBoxJMolScriptString(" boundbox on;");
		          
		      } else if(e.getStateChange()==ItemEvent.DESELECTED){
		          
		    	  theView.setBoxJMolScriptString(" boundbox off;");
		      }
        	
			return;
			
        }  else if( tBtn.equals(theView.getTglBtnPerspective()) ) {        	

		      if(e.getStateChange()==ItemEvent.SELECTED){
		          
		    	  theView.setPerspectiveJMolScriptString(" set perspectiveDepth on;");
		    	  
		      } else if(e.getStateChange()==ItemEvent.DESELECTED){
		    	  
		    	  theView.setPerspectiveJMolScriptString(" set perspectiveDepth off;");
		      }
        	
			return;			
        }
		
	}


	

}

