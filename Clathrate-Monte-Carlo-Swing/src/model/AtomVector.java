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

public class AtomVector {
	
	public double x;
	public double y;
	public double z;
	
	public AtomVector(){
		this.x = 0.0;
		this.y = 0.0;
		this.z = 0.0;
	}
	
	public AtomVector(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public AtomVector(AtomVector v){
		this.x = v.x;
		this.y = v.y;
		this.z = v.z;
	}	
	
	public AtomVector AtomVectorCopy(){
		return (new AtomVector(x, y, z));
	}
	
	public AtomVector AtomVectorAdd(AtomVector v){
		return (new AtomVector(x + v.x, y + v.y, z + v.z));
	}
	
	public AtomVector AtomVectorSubtract(AtomVector v){
		return (new AtomVector(x - v.x, y - v.y, z - v.z));
	}
	
	public AtomVector AtomVectorMultiplyByScalar(double s){
		return (new AtomVector(x * s, y * s, z * s));
	}
	
	public AtomVector AtomVectorDevideByScalar(double s){
		return (new AtomVector(x / s, y / s, z / s));
	}

	public double AtomVectorLength(){
		return Math.sqrt(x*x + y*y + z*z);
	}
	
	public double AtomVectorLengthSquared(){
		return (x*x + y*y + z*z);
	}
	
	public double dotProduct(AtomVector v){
		double sum = this.x * v.x + this.y * v.y + this.z * v.z;
		return sum;
	}
	
	public AtomVector crossProduct(AtomVector v){
		AtomVector prod = new AtomVector(this.y * v.z - this.z * v.y,
										 this.z * v.x - this.x * v.z,
										 this.x * v.y - this.y * v.x
										);
		return prod;
	}
	
	public String toString()
	{
		return "(" + x + ", " + y + ", " + z + ")";
	}
	
	public String getStringXYZ()
	{
		return x + "\t" + y + "\t" + z;
	}
	
	public double getX()
	{
		return x;
	}
	
	public double getY()
	{
		return y;
	}
	
	public double getZ()
	{
		return z;
	}
	
}
