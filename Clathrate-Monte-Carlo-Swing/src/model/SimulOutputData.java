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

import java.awt.geom.Point2D;
import java.util.ArrayList;

public class SimulOutputData {
	
	public int i;
	public double energy;
	public AtomVector box;
	public double energyH;
	public ArrayList<Point2D> gr;
	public StringBuilder sbAllAtomsXYZ;
	
	public void setData(int i, double energy, double energyH, AtomVector box,  ArrayList<Point2D> gr, StringBuilder sbAllAtomsXYZ){
		
		this.i = i;
		this.energy = energy;
		this.box = box;
		this.energyH = energyH;
		this.gr = gr;
		this.sbAllAtomsXYZ = sbAllAtomsXYZ;
	}
	
}
