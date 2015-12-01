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

package view;

import java.awt.EventQueue;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import javax.swing.JDialog;
import javax.swing.BoxLayout;
import javax.swing.GroupLayout;
import javax.swing.GroupLayout.Alignment;
import javax.swing.JPanel;
import javax.swing.LayoutStyle.ComponentPlacement;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;

public class GraphsDialog extends JDialog {
	
	
	private static ArrayList<Point2D> rdf = new ArrayList<Point2D>();
	private static String dialogTitle, yAxesTitle;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					GraphsDialog dialog = new GraphsDialog(rdf, dialogTitle, yAxesTitle);
					dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
					dialog.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the dialog.
	 */
	public GraphsDialog(ArrayList<Point2D> gr, String dialogTitle, String yAxesTitle) {
		
		this.rdf = gr;
		setTitle(dialogTitle);
		setModal(true);
		setBounds(100, 100, 656, 502);
		
		JPanel panelGraphBottom = new JPanel();
		GroupLayout groupLayout = new GroupLayout(getContentPane());
		groupLayout.setHorizontalGroup(
			groupLayout.createParallelGroup(Alignment.TRAILING)
				.addComponent(panelGraphBottom, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 648, Short.MAX_VALUE)
		);
		groupLayout.setVerticalGroup(
			groupLayout.createParallelGroup(Alignment.LEADING)
				.addComponent(panelGraphBottom, GroupLayout.DEFAULT_SIZE, 474, Short.MAX_VALUE)
		);
		panelGraphBottom.setLayout(new BorderLayout(0, 0));
		getContentPane().setLayout(groupLayout);
		
 	
        panelGraphBottom.removeAll();
        panelGraphBottom.revalidate(); 
   	    
 		JFreeChart chartBottom;    
 		XYPlot plotBottom;
 		
		XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries series = new XYSeries("Temp");
        
        dataset.removeAllSeries();
        
        for (int j = 0; j < gr.size(); j++) {
        	series.add(gr.get(j).getX(), gr.get(j).getY());
		}        
        dataset.addSeries(series);
 		
 	    chartBottom = ChartFactory.createXYLineChart(
 	               "", "r, A", yAxesTitle, dataset,
 	               PlotOrientation.VERTICAL, false, false, false);
 	    chartBottom.setBackgroundPaint(Color.white);
 	    
 	    plotBottom = chartBottom.getXYPlot();
 	    
 	    plotBottom.setBackgroundPaint(Color.lightGray);
// 	    plotBottom.setBackgroundPaint(Color.white);
 	    plotBottom.setDomainGridlinePaint(Color.white);
 	    plotBottom.setRangeGridlinePaint(Color.white);
 	    
 	    ValueAxis axis = plotBottom.getDomainAxis();
 	    axis = plotBottom.getRangeAxis();
// 	    axis.setRange(0.0, 4.5);
 	    
 	    XYItemRenderer renderer = plotBottom.getRenderer();
 	    renderer.setSeriesPaint(0, Color.BLUE);
 	    
    
 	    ChartPanel chartPanelBottom;   
 	    panelGraphBottom.setLayout(new BorderLayout(0, 0));
 	    
 	    chartPanelBottom = new ChartPanel(chartBottom) {
 	            @Override
 	            public Dimension getPreferredSize() {
 	                return new Dimension(100, 100);
 	            }
 	    };

 	   chartPanelBottom.setMouseZoomable(true, true);

 	    panelGraphBottom.setLayout(new BorderLayout()); 
 		panelGraphBottom.add(chartPanelBottom);  
 		chartPanelBottom.setLayout(new BorderLayout(0, 0)); 		 
 		panelGraphBottom.repaint(); // This method makes the new chart appear
 		panelGraphBottom.revalidate(); // This removes the old chart
   	      

	}
}
