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

import java.awt.CardLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Locale;

import javax.swing.ButtonGroup;
import javax.swing.GroupLayout;
import javax.swing.ImageIcon;
import javax.swing.GroupLayout.Alignment;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.LayoutStyle.ComponentPlacement;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;
import javax.swing.border.EtchedBorder;
import javax.swing.border.LineBorder;
import javax.swing.border.TitledBorder;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;

import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolViewer;
import org.jmol.util.Logger;

import info.monitorenter.gui.chart.Chart2D;
import info.monitorenter.gui.chart.IAxis;
import info.monitorenter.gui.chart.IAxisScalePolicy;
import info.monitorenter.gui.chart.ITrace2D;
import info.monitorenter.gui.chart.axis.AAxis;
import info.monitorenter.gui.chart.axis.AxisLinear;
import info.monitorenter.gui.chart.axis.scalepolicy.AxisScalePolicyManualTicks;
import info.monitorenter.gui.chart.labelformatters.LabelFormatterNumber;
import info.monitorenter.gui.chart.rangepolicies.RangePolicyFixedViewport;
import info.monitorenter.gui.chart.traces.Trace2DLtd;
import info.monitorenter.gui.chart.traces.Trace2DReplacing;
import info.monitorenter.gui.chart.traces.Trace2DSimple;
import info.monitorenter.util.Range;


import java.awt.BorderLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JToggleButton;
import java.awt.Toolkit;

public class View extends JFrame {

	private JPanel contentPane;
	private final ButtonGroup buttonGroupEnsemble = new ButtonGroup();
	private JTable tableOutput;
	private JRadioButton rdbtnNVT;
	private JRadioButton rdbtnNPT;
	private JSpinner spinnerTemperatureNPT;
	private JSpinner spinnerPressureNPT;
	private JSpinner spinnerFrequencyNPT;
	private JSpinner spinnerNumMCSteps;
	private JLabel lblSimulProgress;
	private JProgressBar progressBarSimulation;
	private JButton btnRun;
	private JButton btnStop;
	private JButton btnRDF;
	private JButton btnNormedRDF;
	private JMenuItem mntmQuit;
	private JMenuItem mntmAbout;
	private JPanel panelTempPressFreqNPT;
	private JPanel panelTempFreqNVT;
	private JPanel panelCardsEnsembles;
	private int rowIndex = 0;
	private JSpinner spinnerTemperatureNVT;
	private JSpinner spinnerFrequencyNVT;
	
	private Chart2D  chartTop;
	private ITrace2D traceTop;
	
	private Chart2D  chartMiddle;
	private ITrace2D traceMiddle;
	
	private Chart2D  chartBottom;
	private ITrace2D traceBottom;
	
	private JToggleButton tglbtnPerspective;
	private JToggleButton tglbtnBox;
	
	private JPanel panelJMol;    
	private JmolPanel jmolPanel = null;

	private String mainJMolScriptString = "select molecule <= 272; wireframe on; spacefill off;";
	private String boxJMolScriptString = " boundbox off;";
	private String perspectiveJMolScriptString = " set perspectiveDepth on;";		
			 
	
	public View() {
		setIconImage(Toolkit.getDefaultToolkit().getImage(View.class.getResource("/resources/ball1.png")));

		setTitle("Monte Carlo Simulation of sH Clathrate Hydrate by Alexander Atamas");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 1170, 827);
		
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		JMenu mnFile = new JMenu("File");
		mnFile.setMnemonic('F');
		menuBar.add(mnFile);
		
		mntmQuit = new JMenuItem("Quit");
		mntmQuit.setIcon(new ImageIcon(View.class.getResource("/resources/exit_16.png")));
		mntmQuit.setMnemonic('Q');
		mntmQuit.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_Q, InputEvent.CTRL_MASK));
		mnFile.add(mntmQuit);
		
		JMenu mnHelp = new JMenu("Help");
		mnHelp.setMnemonic('H');
		menuBar.add(mnHelp);
		
		mntmAbout = new JMenuItem("About");
		mntmAbout.setMnemonic('A');
		mntmAbout.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_MASK));
		mnHelp.add(mntmAbout);
		contentPane = new JPanel();
		contentPane.setAlignmentX(Component.RIGHT_ALIGNMENT);
		contentPane.setAlignmentY(Component.TOP_ALIGNMENT);
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		
		JSplitPane splitMainHorizPane = new JSplitPane();
		splitMainHorizPane.setDividerSize(6);
		splitMainHorizPane.setDividerLocation(300);
		splitMainHorizPane.setAlignmentY(Component.CENTER_ALIGNMENT);
		splitMainHorizPane.setAlignmentX(Component.CENTER_ALIGNMENT);
		
		JPanel panelLeftControlPanel = new JPanel();
		splitMainHorizPane.setLeftComponent(panelLeftControlPanel);
		
		JPanel panelEnsemble = new JPanel();
		panelEnsemble.setBorder(new TitledBorder(new LineBorder(new Color(184, 207, 229)), "Ensemble", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(51, 51, 51)));
		
		JPanel panelMCSteps = new JPanel();
		
		JSeparator separator = new JSeparator();
		
		JPanel panelSimRunStop = new JPanel();
		panelSimRunStop.setBorder(new TitledBorder(new LineBorder(new Color(184, 207, 229)), "Simulation", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(51, 51, 51)));
		
		panelCardsEnsembles = new JPanel();
		GroupLayout gl_panelLeftControlPanel = new GroupLayout(panelLeftControlPanel);
		gl_panelLeftControlPanel.setHorizontalGroup(
			gl_panelLeftControlPanel.createParallelGroup(Alignment.LEADING)
				.addComponent(panelEnsemble, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
				.addComponent(panelCardsEnsembles, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
				.addComponent(panelMCSteps, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
				.addComponent(panelSimRunStop, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
				.addComponent(separator, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
		);
		gl_panelLeftControlPanel.setVerticalGroup(
			gl_panelLeftControlPanel.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelLeftControlPanel.createSequentialGroup()
					.addComponent(panelEnsemble, GroupLayout.PREFERRED_SIZE, 47, GroupLayout.PREFERRED_SIZE)
					.addGap(1)
					.addComponent(panelCardsEnsembles, GroupLayout.PREFERRED_SIZE, 82, GroupLayout.PREFERRED_SIZE)
					.addGap(1)
					.addComponent(separator, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(panelMCSteps, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(panelSimRunStop, GroupLayout.PREFERRED_SIZE, 111, GroupLayout.PREFERRED_SIZE)
					.addContainerGap(462, Short.MAX_VALUE))
		);
		panelCardsEnsembles.setLayout(new CardLayout(0, 0));
		
		panelTempPressFreqNPT = new JPanel();		
		panelCardsEnsembles.add(panelTempPressFreqNPT, "namePanelTempPressFreqNPT");
		panelTempPressFreqNPT.setVisible(false);
		
		JPanel panelTempNPT = new JPanel();
		
		JPanel panelPressureNPT = new JPanel();
		
		JPanel panelFrequencyNPT = new JPanel();
		panelTempPressFreqNPT.setLayout(new GridLayout(0, 1, 0, 0));
		panelTempPressFreqNPT.add(panelTempNPT);
		
		JLabel lblTemperatureK = new JLabel("Temperature, K");
		
		spinnerTemperatureNPT = new JSpinner();
		spinnerTemperatureNPT.setModel(new SpinnerNumberModel(200.0, 10.0, 400.0, 10.0));
		GroupLayout gl_panelTempNPT = new GroupLayout(panelTempNPT);
		gl_panelTempNPT.setHorizontalGroup(
			gl_panelTempNPT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelTempNPT.createSequentialGroup()
					.addGap(12)
					.addComponent(lblTemperatureK, GroupLayout.PREFERRED_SIZE, 125, GroupLayout.PREFERRED_SIZE)
					.addGap(36)
					.addComponent(spinnerTemperatureNPT, GroupLayout.DEFAULT_SIZE, 108, Short.MAX_VALUE)
					.addGap(21))
		);
		gl_panelTempNPT.setVerticalGroup(
			gl_panelTempNPT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelTempNPT.createSequentialGroup()
					.addGroup(gl_panelTempNPT.createParallelGroup(Alignment.LEADING)
						.addGroup(gl_panelTempNPT.createSequentialGroup()
							.addGap(2)
							.addComponent(lblTemperatureK))
						.addComponent(spinnerTemperatureNPT, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		panelTempNPT.setLayout(gl_panelTempNPT);
		panelTempPressFreqNPT.add(panelPressureNPT);
		
		JLabel lblPressureKbars = new JLabel("Pressure, kbars");
		
		spinnerPressureNPT = new JSpinner();
		spinnerPressureNPT.setModel(new SpinnerNumberModel(1.0, 1.0, 50.0, 1.0));
		GroupLayout gl_panelPressureNPT = new GroupLayout(panelPressureNPT);
		gl_panelPressureNPT.setHorizontalGroup(
			gl_panelPressureNPT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelPressureNPT.createSequentialGroup()
					.addGap(12)
					.addComponent(lblPressureKbars, GroupLayout.PREFERRED_SIZE, 149, GroupLayout.PREFERRED_SIZE)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(spinnerPressureNPT, GroupLayout.DEFAULT_SIZE, 108, Short.MAX_VALUE)
					.addGap(21))
		);
		gl_panelPressureNPT.setVerticalGroup(
			gl_panelPressureNPT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelPressureNPT.createSequentialGroup()
					.addGroup(gl_panelPressureNPT.createParallelGroup(Alignment.LEADING)
						.addGroup(gl_panelPressureNPT.createSequentialGroup()
							.addGap(3)
							.addComponent(lblPressureKbars))
						.addComponent(spinnerPressureNPT, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		panelPressureNPT.setLayout(gl_panelPressureNPT);
		panelTempPressFreqNPT.add(panelFrequencyNPT);
		
		JLabel labelFrequency = new JLabel("Frequency of output");
		
		spinnerFrequencyNPT = new JSpinner();
		spinnerFrequencyNPT.setModel(new SpinnerNumberModel(2000, 500, 100000, 200));
		GroupLayout gl_panelFrequencyNPT = new GroupLayout(panelFrequencyNPT);
		gl_panelFrequencyNPT.setHorizontalGroup(
			gl_panelFrequencyNPT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelFrequencyNPT.createSequentialGroup()
					.addGap(12)
					.addComponent(labelFrequency, GroupLayout.PREFERRED_SIZE, 149, GroupLayout.PREFERRED_SIZE)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(spinnerFrequencyNPT, GroupLayout.DEFAULT_SIZE, 108, Short.MAX_VALUE)
					.addGap(21))
		);
		gl_panelFrequencyNPT.setVerticalGroup(
			gl_panelFrequencyNPT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelFrequencyNPT.createSequentialGroup()
					.addGroup(gl_panelFrequencyNPT.createParallelGroup(Alignment.LEADING)
						.addGroup(gl_panelFrequencyNPT.createSequentialGroup()
							.addGap(2)
							.addComponent(labelFrequency))
						.addComponent(spinnerFrequencyNPT, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		panelFrequencyNPT.setLayout(gl_panelFrequencyNPT);
		
		panelTempFreqNVT = new JPanel();
		panelCardsEnsembles.add(panelTempFreqNVT, "namePanelTempFreqNVT");
		panelTempFreqNVT.setVisible(true);
		
		JPanel panelTempNVT = new JPanel();
		
		JLabel label_1 = new JLabel("Temperature, K");
		
		spinnerTemperatureNVT = new JSpinner();
		spinnerTemperatureNVT.setModel(new SpinnerNumberModel(200.0, 10.0, 400.0, 10.0));
		GroupLayout gl_panelTempNVT = new GroupLayout(panelTempNVT);
		gl_panelTempNVT.setHorizontalGroup(
			gl_panelTempNVT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelTempNVT.createSequentialGroup()
					.addContainerGap()
					.addComponent(label_1, GroupLayout.PREFERRED_SIZE, 125, GroupLayout.PREFERRED_SIZE)
					.addGap(36)
					.addComponent(spinnerTemperatureNVT, GroupLayout.DEFAULT_SIZE, 91, Short.MAX_VALUE)
					.addGap(21))
		);
		gl_panelTempNVT.setVerticalGroup(
			gl_panelTempNVT.createParallelGroup(Alignment.LEADING)
				.addGroup(Alignment.TRAILING, gl_panelTempNVT.createSequentialGroup()
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
					.addGroup(gl_panelTempNVT.createParallelGroup(Alignment.LEADING)
						.addGroup(gl_panelTempNVT.createSequentialGroup()
							.addGap(2)
							.addComponent(label_1))
						.addComponent(spinnerTemperatureNVT, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
					.addContainerGap())
		);
		panelTempNVT.setLayout(gl_panelTempNVT);
		
		JPanel panelFreqNVT = new JPanel();
		
		JLabel label_3 = new JLabel("Frequency of output");
		
		spinnerFrequencyNVT = new JSpinner();
		spinnerFrequencyNVT.setModel(new SpinnerNumberModel(3000, 500, 100000, 200));
		GroupLayout gl_panelFreqNVT = new GroupLayout(panelFreqNVT);
		gl_panelFreqNVT.setHorizontalGroup(
			gl_panelFreqNVT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelFreqNVT.createSequentialGroup()
					.addContainerGap()
					.addComponent(label_3, GroupLayout.PREFERRED_SIZE, 149, GroupLayout.PREFERRED_SIZE)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(spinnerFrequencyNVT, GroupLayout.DEFAULT_SIZE, 91, Short.MAX_VALUE)
					.addGap(21))
		);
		gl_panelFreqNVT.setVerticalGroup(
			gl_panelFreqNVT.createParallelGroup(Alignment.LEADING)
				.addGroup(Alignment.TRAILING, gl_panelFreqNVT.createSequentialGroup()
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
					.addGroup(gl_panelFreqNVT.createParallelGroup(Alignment.LEADING)
						.addGroup(gl_panelFreqNVT.createSequentialGroup()
							.addGap(2)
							.addComponent(label_3))
						.addComponent(spinnerFrequencyNVT, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
					.addContainerGap())
		);
		panelFreqNVT.setLayout(gl_panelFreqNVT);
		GroupLayout gl_panelTempFreqNVT = new GroupLayout(panelTempFreqNVT);
		gl_panelTempFreqNVT.setHorizontalGroup(
			gl_panelTempFreqNVT.createParallelGroup(Alignment.LEADING)
				.addComponent(panelTempNVT, GroupLayout.DEFAULT_SIZE, 336, Short.MAX_VALUE)
				.addComponent(panelFreqNVT, GroupLayout.DEFAULT_SIZE, 336, Short.MAX_VALUE)
		);
		gl_panelTempFreqNVT.setVerticalGroup(
			gl_panelTempFreqNVT.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelTempFreqNVT.createSequentialGroup()
					.addComponent(panelTempNVT, GroupLayout.PREFERRED_SIZE, 41, GroupLayout.PREFERRED_SIZE)
					.addGap(1)
					.addComponent(panelFreqNVT, GroupLayout.PREFERRED_SIZE, 41, GroupLayout.PREFERRED_SIZE)
					.addContainerGap())
		);
		panelTempFreqNVT.setLayout(gl_panelTempFreqNVT);
		
		progressBarSimulation = new JProgressBar(0, 100);
		progressBarSimulation.setValue(0);
		progressBarSimulation.setStringPainted(true);
		
		lblSimulProgress = new JLabel("Progress:");
		
		btnRun = new JButton("Run");
		
		btnStop = new JButton("Stop");
		btnStop.setEnabled(false);
		GroupLayout gl_panelSimRunStop = new GroupLayout(panelSimRunStop);
		gl_panelSimRunStop.setHorizontalGroup(
			gl_panelSimRunStop.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelSimRunStop.createSequentialGroup()
					.addContainerGap()
					.addGroup(gl_panelSimRunStop.createParallelGroup(Alignment.LEADING)
						.addComponent(progressBarSimulation, GroupLayout.DEFAULT_SIZE, 251, Short.MAX_VALUE)
						.addComponent(lblSimulProgress)
						.addGroup(gl_panelSimRunStop.createSequentialGroup()
							.addComponent(btnRun)
							.addPreferredGap(ComponentPlacement.RELATED)
							.addComponent(btnStop)))
					.addContainerGap())
		);
		gl_panelSimRunStop.setVerticalGroup(
			gl_panelSimRunStop.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelSimRunStop.createSequentialGroup()
					.addGap(6)
					.addComponent(lblSimulProgress)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(progressBarSimulation, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addGroup(gl_panelSimRunStop.createParallelGroup(Alignment.BASELINE)
						.addComponent(btnRun)
						.addComponent(btnStop))
					.addContainerGap(17, Short.MAX_VALUE))
		);
		panelSimRunStop.setLayout(gl_panelSimRunStop);
		
		JLabel lblNumberMCSteps = new JLabel("Number of MC steps");
		
		spinnerNumMCSteps = new JSpinner();
		spinnerNumMCSteps.setModel(new SpinnerNumberModel(250000, 10000, 10000000, 25000));
		GroupLayout gl_panelMCSteps = new GroupLayout(panelMCSteps);
		gl_panelMCSteps.setHorizontalGroup(
			gl_panelMCSteps.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelMCSteps.createSequentialGroup()
					.addContainerGap()
					.addComponent(lblNumberMCSteps)
					.addPreferredGap(ComponentPlacement.UNRELATED)
					.addComponent(spinnerNumMCSteps, GroupLayout.PREFERRED_SIZE, 106, Short.MAX_VALUE)
					.addGap(21))
		);
		gl_panelMCSteps.setVerticalGroup(
			gl_panelMCSteps.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelMCSteps.createSequentialGroup()
					.addGap(9)
					.addGroup(gl_panelMCSteps.createParallelGroup(Alignment.BASELINE)
						.addComponent(lblNumberMCSteps)
						.addComponent(spinnerNumMCSteps, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
					.addContainerGap(14, Short.MAX_VALUE))
		);
		panelMCSteps.setLayout(gl_panelMCSteps);
		
		rdbtnNVT = new JRadioButton("NVT");
		buttonGroupEnsemble.add(rdbtnNVT);
		
		rdbtnNPT = new JRadioButton("NPT");
		rdbtnNPT.setSelected(true);
		buttonGroupEnsemble.add(rdbtnNPT);
		GroupLayout gl_panelEnsemble = new GroupLayout(panelEnsemble);
		gl_panelEnsemble.setHorizontalGroup(
			gl_panelEnsemble.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelEnsemble.createSequentialGroup()
					.addComponent(rdbtnNVT)
					.addGap(18)
					.addComponent(rdbtnNPT)
					.addContainerGap(204, Short.MAX_VALUE))
		);
		gl_panelEnsemble.setVerticalGroup(
			gl_panelEnsemble.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelEnsemble.createSequentialGroup()
					.addGroup(gl_panelEnsemble.createParallelGroup(Alignment.BASELINE)
						.addComponent(rdbtnNVT)
						.addComponent(rdbtnNPT))
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		panelEnsemble.setLayout(gl_panelEnsemble);
		panelLeftControlPanel.setLayout(gl_panelLeftControlPanel);
		GroupLayout gl_contentPane = new GroupLayout(contentPane);
		gl_contentPane.setHorizontalGroup(
			gl_contentPane.createParallelGroup(Alignment.LEADING)
				.addComponent(splitMainHorizPane, GroupLayout.DEFAULT_SIZE, 880, Short.MAX_VALUE)
		);
		gl_contentPane.setVerticalGroup(
			gl_contentPane.createParallelGroup(Alignment.LEADING)
				.addComponent(splitMainHorizPane, Alignment.TRAILING, GroupLayout.DEFAULT_SIZE, 678, Short.MAX_VALUE)
		);
		
		JPanel panelRightMainPanel = new JPanel();
		splitMainHorizPane.setRightComponent(panelRightMainPanel);
		
		JSplitPane splitPaneVertOutputPanel = new JSplitPane();
		splitPaneVertOutputPanel.setDividerSize(6);
		splitPaneVertOutputPanel.setResizeWeight(1);
		splitPaneVertOutputPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);
		
		splitPaneVertOutputPanel.setDividerLocation(605);
		splitPaneVertOutputPanel.setMinimumSize( new Dimension( 250, 400 ) );
		splitPaneVertOutputPanel.setPreferredSize( new Dimension( 250, 200 ) );
		
		GroupLayout gl_panelRightMainPanel = new GroupLayout(panelRightMainPanel);
		gl_panelRightMainPanel.setHorizontalGroup(
			gl_panelRightMainPanel.createParallelGroup(Alignment.LEADING)
				.addComponent(splitPaneVertOutputPanel, Alignment.TRAILING, GroupLayout.DEFAULT_SIZE, 571, Short.MAX_VALUE)
		);
		gl_panelRightMainPanel.setVerticalGroup(
			gl_panelRightMainPanel.createParallelGroup(Alignment.LEADING)
				.addComponent(splitPaneVertOutputPanel, GroupLayout.DEFAULT_SIZE, 676, Short.MAX_VALUE)
		);
		
		JPanel panelGraphicalOutput = new JPanel();
		panelGraphicalOutput.setPreferredSize(new Dimension(10, 100));
		splitPaneVertOutputPanel.setLeftComponent(panelGraphicalOutput);
		
		JSplitPane splitPaneHorizGraphicalOutput = new JSplitPane();
		splitPaneHorizGraphicalOutput.setDividerSize(6);
		splitPaneHorizGraphicalOutput.setDividerLocation(500);
		splitPaneHorizGraphicalOutput.setResizeWeight(1); 
		
		
		GroupLayout gl_panelGraphicalOutput = new GroupLayout(panelGraphicalOutput);
		gl_panelGraphicalOutput.setHorizontalGroup(
			gl_panelGraphicalOutput.createParallelGroup(Alignment.LEADING)
				.addComponent(splitPaneHorizGraphicalOutput, GroupLayout.DEFAULT_SIZE, 585, Short.MAX_VALUE)
		);
		gl_panelGraphicalOutput.setVerticalGroup(
			gl_panelGraphicalOutput.createParallelGroup(Alignment.LEADING)
				.addComponent(splitPaneHorizGraphicalOutput, GroupLayout.DEFAULT_SIZE, 499, Short.MAX_VALUE)
		);
		
		JPanel panelMolViewer = new JPanel();
		splitPaneHorizGraphicalOutput.setLeftComponent(panelMolViewer);
		
		JPanel panelJMolToolbar = new JPanel();
		panelJMolToolbar.setBorder(new LineBorder(new Color(0, 0, 0)));
		
		panelJMol = new JPanel();
		
		GroupLayout gl_panelMolViewer = new GroupLayout(panelMolViewer);
		gl_panelMolViewer.setHorizontalGroup(
			gl_panelMolViewer.createParallelGroup(Alignment.LEADING)
				.addComponent(panelJMolToolbar, GroupLayout.DEFAULT_SIZE, 559, Short.MAX_VALUE)
				.addComponent(panelJMol, GroupLayout.DEFAULT_SIZE, 559, Short.MAX_VALUE)
		);
		gl_panelMolViewer.setVerticalGroup(
			gl_panelMolViewer.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelMolViewer.createSequentialGroup()
					.addComponent(panelJMolToolbar, GroupLayout.PREFERRED_SIZE, 44, GroupLayout.PREFERRED_SIZE)
					.addGap(1)
					.addComponent(panelJMol, GroupLayout.DEFAULT_SIZE, 490, Short.MAX_VALUE))
		);
		
		tglbtnBox = new JToggleButton("Show Box");
		tglbtnBox.setToolTipText("Turn on/off bounding box");
		
		tglbtnPerspective = new JToggleButton("Show Perspective");
		tglbtnPerspective.setSelected(true);
		tglbtnPerspective.setToolTipText("Turn on/off perspective depth");
		GroupLayout gl_panelJMolToolbar = new GroupLayout(panelJMolToolbar);
		gl_panelJMolToolbar.setHorizontalGroup(
			gl_panelJMolToolbar.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelJMolToolbar.createSequentialGroup()
					.addContainerGap()
					.addComponent(tglbtnBox)
					.addPreferredGap(ComponentPlacement.RELATED)
					.addComponent(tglbtnPerspective)
					.addContainerGap(305, Short.MAX_VALUE))
		);
		gl_panelJMolToolbar.setVerticalGroup(
			gl_panelJMolToolbar.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelJMolToolbar.createSequentialGroup()
					.addGap(10)
					.addGroup(gl_panelJMolToolbar.createParallelGroup(Alignment.BASELINE)
						.addComponent(tglbtnBox)
						.addComponent(tglbtnPerspective))
					.addContainerGap(GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
		);
		panelJMolToolbar.setLayout(gl_panelJMolToolbar);
		panelJMol.setLayout(new BorderLayout(0, 0));
		
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		jmolPanel = new JmolPanel();
		jmolPanel.setPreferredSize(new Dimension(400, 400));
		panelJMol.add(jmolPanel);
/*		String strError;
		strError = jmolPanel.viewer.openStringInline(strXyzHOH);*/
		
	    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		panelMolViewer.setLayout(gl_panelMolViewer);
		
		JPanel panelGraphOutput = new JPanel();
		splitPaneHorizGraphicalOutput.setRightComponent(panelGraphOutput);
		
		JPanel panelGraphTop = new JPanel();
		panelGraphTop.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		
		JPanel panelGraphMiddle = new JPanel();
		panelGraphMiddle.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		
		JPanel panelGraphBottom = new JPanel();
		panelGraphBottom.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		panelGraphTop.setLayout(new BorderLayout(0, 0));
		panelGraphBottom.setLayout(new BorderLayout(0, 0));
		panelGraphMiddle.setLayout(new BorderLayout(0, 0));
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		
		Font titleFont;
		IAxis.AxisTitle axisTitleY;
		IAxis.AxisTitle axisTitleX;
		AAxis<IAxisScalePolicy> axisXTop;
		AAxis<IAxisScalePolicy> axisYRight;
		
		
		chartTop = new Chart2D();
//		traceTop = new Trace2DLtd(100);
		traceTop = new Trace2DSimple();
		traceTop.setColor(Color.BLUE);
		traceTop.setName("");
		chartTop.addTrace(traceTop);
	
		
		IAxis<IAxisScalePolicy> yAxis = (IAxis<IAxisScalePolicy>)this.chartTop.getAxisY();
		yAxis.setRangePolicy(new RangePolicyFixedViewport(new Range(-10, -21)));
		yAxis.setAxisScalePolicy(new AxisScalePolicyManualTicks());
		yAxis.setMinorTickSpacing(3.0);
		yAxis.setFormatter(new LabelFormatterNumber(new DecimalFormat("0.00")));
		
	    titleFont = UIManager.getDefaults().getFont("Label.font").deriveFont(11f).deriveFont(Font.BOLD);
		
	    axisTitleY = chartTop.getAxisY().getAxisTitle();
	    axisTitleY.setTitle("U, kcal/mol");
	    axisTitleY.setTitleFont(titleFont);
	    
	    axisTitleX = chartTop.getAxisX().getAxisTitle();
	    axisTitleX.setTitle("MC step");
	    axisTitleX.setTitleFont(titleFont);
	    
	    axisXTop = new AxisLinear<IAxisScalePolicy>();
	    axisXTop.setPaintScale(false);
	    axisYRight = new AxisLinear<IAxisScalePolicy>();
	    axisYRight.setPaintScale(false);
	    chartTop.setAxisXTop(axisXTop, 0);
	    chartTop.setAxisYRight(axisYRight,0);		
		chartTop.getAxisX().setPaintGrid(true);
		chartTop.getAxisY().setPaintGrid(true);   
		
		panelGraphTop.setLayout(new BorderLayout(0, 0));
	    panelGraphTop.add(chartTop);
//		panelJMol.add(chartTop);
	    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	    
	    
	    chartMiddle = new Chart2D();
//		traceMiddle = new Trace2DLtd(200);
	    traceMiddle = new Trace2DSimple();
	    traceMiddle.setColor(Color.GREEN);
	    traceMiddle.setName("");
		chartMiddle.addTrace(traceMiddle);
		
		yAxis = (IAxis<IAxisScalePolicy>)this.chartMiddle.getAxisY();
		yAxis.setRangePolicy(new RangePolicyFixedViewport(new Range(0.75, 0.85)));
		yAxis.setAxisScalePolicy(new AxisScalePolicyManualTicks());
		yAxis.setMinorTickSpacing(0.025);
		yAxis.setFormatter(new LabelFormatterNumber(new DecimalFormat("0.00")));
		
	    titleFont = UIManager.getDefaults().getFont("Label.font").deriveFont(11f).deriveFont(Font.BOLD);
		
	    axisTitleY = chartMiddle.getAxisY().getAxisTitle();
	    axisTitleY.setTitle("rho, g/cm^3");
	    axisTitleY.setTitleFont(titleFont);
	    
	    axisTitleX = chartMiddle.getAxisX().getAxisTitle();
	    axisTitleX.setTitle("MC step");
	    axisTitleX.setTitleFont(titleFont);
	    
	    axisXTop = new AxisLinear<IAxisScalePolicy>();
	    axisXTop.setPaintScale(false);
	    axisYRight = new AxisLinear<IAxisScalePolicy>();
	    axisYRight.setPaintScale(false);
	    chartMiddle.setAxisXTop(axisXTop, 0);
	    chartMiddle.setAxisYRight(axisYRight,0);
	    chartMiddle.getAxisX().setPaintGrid(true);
	    chartMiddle.getAxisY().setPaintGrid(true); 
	    
		panelGraphMiddle.setLayout(new BorderLayout(0, 0));
		panelGraphMiddle.add(chartMiddle);			
//		panelJMol.add(chartMiddle);	
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		chartBottom = new Chart2D();
		traceBottom = new Trace2DSimple();
//		traceBottom = new Trace2DLtd();
		traceBottom.setColor(Color.RED);
		traceBottom.setName("");
	    chartBottom.addTrace(traceBottom);
		
		yAxis = (IAxis<IAxisScalePolicy>)this.chartBottom.getAxisY();
		yAxis.setRangePolicy(new RangePolicyFixedViewport(new Range(24.3, 24.76)));
		yAxis.setAxisScalePolicy(new AxisScalePolicyManualTicks());
		yAxis.setMinorTickSpacing(0.15);
		yAxis.setFormatter(new LabelFormatterNumber(new DecimalFormat("0.00")));

		
	    titleFont = UIManager.getDefaults().getFont("Label.font").deriveFont(11f).deriveFont(Font.BOLD);
		
	    axisTitleY = chartBottom.getAxisY().getAxisTitle();
	    axisTitleY.setTitle("box, A");
	    axisTitleY.setTitleFont(titleFont);
	    
	    axisTitleX = chartBottom.getAxisX().getAxisTitle();
	    axisTitleX.setTitle("MC step");
	    axisTitleX.setTitleFont(titleFont);
	    
	    axisXTop = new AxisLinear<IAxisScalePolicy>();
	    axisXTop.setPaintScale(false);
	    axisYRight = new AxisLinear<IAxisScalePolicy>();
	    axisYRight.setPaintScale(false);
	    chartBottom.setAxisXTop(axisXTop, 0);
	    chartBottom.setAxisYRight(axisYRight,0);
	    chartBottom.getAxisX().setPaintGrid(true);
	    chartBottom.getAxisY().setPaintGrid(true); 
	    
		panelGraphBottom.setLayout(new BorderLayout(0, 0));
		panelGraphBottom.add(chartBottom);
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		

		
		JPanel panel = new JPanel();
		panel.setBorder(new LineBorder(new Color(0, 0, 0)));
		GroupLayout gl_panelGraphOutput = new GroupLayout(panelGraphOutput);
		gl_panelGraphOutput.setHorizontalGroup(
			gl_panelGraphOutput.createParallelGroup(Alignment.TRAILING)
				.addGroup(gl_panelGraphOutput.createSequentialGroup()
					.addGroup(gl_panelGraphOutput.createParallelGroup(Alignment.TRAILING)
						.addComponent(panelGraphBottom, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
						.addComponent(panelGraphTop, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
						.addComponent(panel, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
						.addComponent(panelGraphMiddle, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE))
					.addGap(0))
		);
		gl_panelGraphOutput.setVerticalGroup(
			gl_panelGraphOutput.createParallelGroup(Alignment.LEADING)
				.addGroup(gl_panelGraphOutput.createSequentialGroup()
					.addComponent(panel, GroupLayout.PREFERRED_SIZE, 43, GroupLayout.PREFERRED_SIZE)
					.addGap(2)
					.addComponent(panelGraphTop, GroupLayout.DEFAULT_SIZE, 187, Short.MAX_VALUE)
					.addGap(2)
					.addComponent(panelGraphMiddle, GroupLayout.DEFAULT_SIZE, 179, Short.MAX_VALUE)
					.addGap(2)
					.addComponent(panelGraphBottom, GroupLayout.DEFAULT_SIZE, 187, Short.MAX_VALUE))
		);
		
		btnRDF = new JButton("g(r)");
		btnRDF.setToolTipText("Shows \"Radial Distribution Function (RDF)\" graph");
		btnRDF.setEnabled(false);
		btnRDF.setBounds(13, 9, 59, 25);
		panel.setLayout(null);
		
		btnNormedRDF = new JButton("gn(r)");
		btnNormedRDF.setToolTipText("Shows \"Particle distribution\" graph ");
		btnNormedRDF.setEnabled(false);
		btnNormedRDF.setBounds(78, 9, 68, 25);
		panel.add(btnNormedRDF);
		panel.add(btnRDF);
		panelGraphOutput.setLayout(gl_panelGraphOutput);
		panelGraphicalOutput.setLayout(gl_panelGraphicalOutput);
		
		JScrollPane scrollPaneTable = new JScrollPane();
		splitPaneVertOutputPanel.setRightComponent(scrollPaneTable);
		
		tableOutput = new JTable();
		tableOutput.setModel(new DefaultTableModel(
			new Object[][] {
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
				{null, null, null, null, null, null},
			},
//				new Object[][] {},
			new String[] {
				"MC step", "U, kcal/mol", "H, kcal/mol", "xBox, A", "yBox, A", "zBox, A"
			}
		) {
			Class[] columnTypes = new Class[] {
				String.class, String.class, String.class, String.class, String.class, String.class
			};
			public Class getColumnClass(int columnIndex) {
				return columnTypes[columnIndex];
			}
			boolean[] columnEditables = new boolean[] {
				false, false, false, false, false, false
			};
			public boolean isCellEditable(int row, int column) {
				return columnEditables[column];
			}
		});
		tableOutput.getColumnModel().getColumn(0).setPreferredWidth(55);
		tableOutput.getColumnModel().getColumn(0).setMinWidth(55);
		tableOutput.getColumnModel().getColumn(1).setPreferredWidth(95);
		tableOutput.getColumnModel().getColumn(1).setMinWidth(95);
		tableOutput.getColumnModel().getColumn(2).setPreferredWidth(95);
		tableOutput.getColumnModel().getColumn(2).setMinWidth(95);
		tableOutput.getColumnModel().getColumn(3).setMinWidth(75);
		tableOutput.getColumnModel().getColumn(4).setMinWidth(75);
		tableOutput.getColumnModel().getColumn(5).setMinWidth(75);
		scrollPaneTable.setViewportView(tableOutput);
		panelRightMainPanel.setLayout(gl_panelRightMainPanel);
		contentPane.setLayout(gl_contentPane);
		
		DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
		centerRenderer.setHorizontalAlignment( SwingConstants.CENTER );
		tableOutput.setDefaultRenderer(String.class, centerRenderer);
		
		JTableHeader header = tableOutput.getTableHeader();
		header.setFont(new Font("Dialog", Font.BOLD, 12));
		
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public void quitApp(){
		
		System.exit(0);
	}
	
	public void showAboutDialog(){
		
		About about = new About();		
		about.setModal(true);
		about.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		about.setLocationRelativeTo(this);
		about.setVisible(true);
	}
	
	public void showGraphsDialog(ArrayList<Point2D> gr, String dialogTitle, String yAxesTitle){
		
		GraphsDialog dialog = new GraphsDialog(gr, dialogTitle, yAxesTitle);
		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		dialog.setLocationRelativeTo(this);
		dialog.setVisible(true);
	}
	
		
	public JToggleButton getTglBtnBox(){
		
		return tglbtnBox;
	}
	
	public JToggleButton getTglBtnPerspective(){
		
		return tglbtnPerspective;
	}

	public JProgressBar getProgressBarSimulation(){
		
		return progressBarSimulation;
	}
	
	public JButton getBtnRDF(){
		return btnRDF;
	}
	
	public JButton getBtnNormedRDF(){
		return btnNormedRDF;
	}

	public JButton getBtnRun(){
		rowIndex = 0;
		return btnRun;
	}
	
	public JButton getBtnStop(){
		rowIndex = 0;
		return btnStop;
	}	
	
	public JMenuItem getMntmQuit(){
		return mntmQuit;
	}	

	public JMenuItem getMntmAbout(){
		return mntmAbout;
	}	

	public JRadioButton getRdbtnNVT(){
		return rdbtnNVT;
	}
	
	public JRadioButton getRdbtnNPT(){
		return rdbtnNPT;
	}

	public JPanel getPanelCardsEnsembles(){
		return panelCardsEnsembles;
	}
	
	public int getNumbMCSteps(){
		
		return (int)spinnerNumMCSteps.getValue();
	}
	
	public int getFrequency(){
		
		int val = 0;
		
		if ( rdbtnNVT.isSelected() ){			
			
			val = (int)spinnerFrequencyNVT.getValue();
			
		} else if ( rdbtnNPT.isSelected() ){
			
			val = (int)spinnerFrequencyNPT.getValue();
		}		
		return val;
	}
	
	public double getTemperature(){
		
		double val = 0.0;
		
		if ( rdbtnNVT.isSelected() ){			
			
			val = (double)spinnerTemperatureNVT.getValue();
			
		} else if ( rdbtnNPT.isSelected() ){
			
			val = (double)spinnerTemperatureNPT.getValue();
		}		
		return val;
	}
	
	public double getPressure(){
		
		double val = 0.0;
		
		if (rdbtnNPT.isSelected()) {

			val = (double) spinnerPressureNPT.getValue();
		}
		return val;
	}
	
	public boolean isNVTSelected(){
		
		return rdbtnNVT.isSelected();
	}
	
	public boolean isNPTSelected(){
		
		return rdbtnNPT.isSelected();
	}
	
	public void setBoxJMolScriptString(String jMolScriptString){
		
		boxJMolScriptString = jMolScriptString;
		
		if ( !btnStop.isEnabled() ){
			
			mainJMolScriptString += boxJMolScriptString;
			mainJMolScriptString += perspectiveJMolScriptString;
			jmolPanel.viewer.evalString(mainJMolScriptString);				
		}
	}
	
	public void setPerspectiveJMolScriptString(String jMolScriptString){
		
		perspectiveJMolScriptString = jMolScriptString;
		
		if ( !btnStop.isEnabled() ){
			
			mainJMolScriptString += boxJMolScriptString;
			mainJMolScriptString += perspectiveJMolScriptString;
			jmolPanel.viewer.evalString(mainJMolScriptString);			
		}
	}
	
	public void showCrystalInJMol(StringBuilder sbAllAtomsXYZ){
		
		if ( jmolPanel != null ){

			String strError = jmolPanel.viewer.openStringInline(sbAllAtomsXYZ.toString());
			
			mainJMolScriptString += boxJMolScriptString;
			mainJMolScriptString += perspectiveJMolScriptString;
			
		    if (strError == null){
		        jmolPanel.viewer.evalString(mainJMolScriptString);
		    } else {
		        Logger.error(strError);
		    }
		    
		    try {
		          Thread.sleep(20);
		    } catch (InterruptedException e) {
		    	  System.out.println("Thread.sleep Exception");
		    }
		    
		

		}
	}
	
	public void drawGraphsTitle(){
		
        traceTop.addPoint(0.0, 0.0);
        traceMiddle.addPoint(0.0, 0.0);
        traceBottom.addPoint(0.0, 0.0);		
	}
	
	public void showTempPressFreqNPT(){
		
		CardLayout cl = (CardLayout) (panelCardsEnsembles.getLayout());
		cl.show(panelCardsEnsembles, "namePanelTempPressFreqNPT"); 
	
	}	
	
	public void showTempFreqNVT(){

		CardLayout cl = (CardLayout) (panelCardsEnsembles.getLayout());
		cl.show(panelCardsEnsembles, "namePanelTempFreqNVT"); 
	}
	
	public void setSimulCancelledText(){
		
		lblSimulProgress.setText("Simulation cancelled.");
	}
	
	public void setSimulCompletedText(){
		
		lblSimulProgress.setText("Simulation completed.");	
		btnRun.setEnabled(true);
		btnStop.setEnabled(false);
	}
	
	public void setSimulProgressText(int i, int numbMC){
		
		Locale.setDefault(Locale.US);
		lblSimulProgress.setText(String.format("Progress: %,d out of %,d", i, numbMC));		
	}	
	
	public void emptyOutputTable(){
		 
		DefaultTableModel tableOutputModel = (DefaultTableModel) tableOutput.getModel();
		tableOutputModel.setRowCount(0);
		tableOutputModel.setRowCount(8);
		
		traceTop.removeAllPoints();
		traceMiddle.removeAllPoints();
		traceBottom.removeAllPoints(); 
	}
	
	public void updateGUI(int i, double enU, double enH, double xBox, double yBox, double zBox, ArrayList<Point2D> gr, double density, StringBuilder sbAllAtomsXYZ){
		
		String str_i = String.format("%7d", i);
		String str_enU = String.format("%7.3f", enU);
		String str_enH = String.format("%7.3f", enH);
		String str_xBox = String.format("%7.3f", xBox);
		String str_yBox = String.format("%7.3f", yBox);
		String str_zBox = String.format("%7.3f", zBox);
		
		DefaultTableModel tableOutputModel = (DefaultTableModel) tableOutput.getModel();	
		
		if ( rowIndex > 7 ){
			
			tableOutputModel.addRow(new Object[]{str_i, str_enU, str_enH, str_xBox, str_yBox, str_zBox});		

		} else {
	        tableOutput.setValueAt(str_i, rowIndex, 0);
	        tableOutput.setValueAt(str_enU, rowIndex, 1);
	        tableOutput.setValueAt(str_enH, rowIndex, 2);
	        tableOutput.setValueAt(str_xBox, rowIndex, 3);
	        tableOutput.setValueAt(str_yBox, rowIndex, 4);
	        tableOutput.setValueAt(str_zBox, rowIndex, 5);				
		}
        rowIndex++;		
        
        showCrystalInJMol(sbAllAtomsXYZ);
        
/*		traceTop.removeAllPoints();
		traceMiddle.removeAllPoints();   */  
        traceTop.addPoint((double)i, enU);        
        traceMiddle.addPoint((double)i, density); 
        traceBottom.addPoint((double)i, xBox);
       
//        traceBottom.removeAllPoints();
/*        gr.add(new Point2D.Double(Double.NaN, Double.NaN));
		for (int j = 0; j < gr.size(); j++) {
			traceBottom.addPoint(gr.get(j).getX(), gr.get(j).getY());
		} */  
		
		int numbOfRows = tableOutputModel.getRowCount()-1;	
		tableOutput.getSelectionModel().setSelectionInterval(numbOfRows, numbOfRows);
		tableOutput.scrollRectToVisible(new Rectangle(tableOutput.getCellRect(numbOfRows, 0, true)));		
	}	
	
	  private final static String strXyzHOH = 
		      "3\n" +
		      "water\n" +
		  		"O  0.0 0.0 0.0\n" +
		  		"H  0.76923955 -0.59357141 0.0\n" +
		  		"H -0.76923955 -0.59357141 0.0\n";

	  
	  static class JmolPanel extends JPanel {

		    JmolViewer viewer;
		    
		    private final Dimension currentSize = new Dimension();
		    
		    JmolPanel() {
		      viewer = JmolViewer.allocateViewer(this, new SmarterJmolAdapter(), 
		          null, null, null, null, null);
		    }

		    @Override
		    public void paint(Graphics g) {
		      getSize(currentSize);
		      viewer.renderScreenImage(g, currentSize.width, currentSize.height);
		    }
		  }
}




