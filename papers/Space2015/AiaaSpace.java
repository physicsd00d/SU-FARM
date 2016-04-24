package aiaaSpace;

/* 
 * This package is used to be a clean run through of the level zero verification.
 * 		All prototyping occurs elsewhere.  There should be nothing extraneous in here.
 */

// ========== FACET imports
import facet_pkg.FACETInterface;
import facet_pkg.RcdAcPosSetup;
import facet_pkg.UserSUAListPnl;
import facet_pkg.api.server.*;	// This is the main way we communicate with facet, through the server
import facet_pkg.api.server_api.*;

import facet_pkg.api.server_api.SUAInterface;
import facet_pkg.avoidance.Avoidance;
import facet_pkg.avoidance.Avoidance.AvoidanceAlgorithm;

import facet_pkg.aircraft_filter.AircraftFilterDialog;
import facet_pkg.aircraft_filter.AircraftFilterManager;
import facet_pkg.aircraft_filter.AircraftFilterPreferences;
import facet_pkg.aircraft_filter.AircraftFilterData;
import facet_pkg.aircraft_filter.FilteredAircraftDisplay;
import facet_pkg.aircraft_filter.FilteredAircraftTableManager;
import facet_pkg.aircraft_filter.FilteredAircraftTableManager.DumpException;

import facet_pkg.aircraft_filter.FilteredAircraftDisplay;

//// ========== JAVA imports
import java.io.File;			// Needed in order to load folders / files
import java.io.FilenameFilter;
import java.util.Arrays;
//import org.apache.commons.io.filefilter.WildcardFileFilter;  //Can't find this!?
import java.io.FileInputStream;
import java.nio.CharBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.util.ArrayList;
import java.util.List;




// For trying to decode ASDI files
//import java.io.FileInputStream;

public class AiaaSpace {
	
	//Initiates connection to FACET
	FACETServerAPI srv = null;
	
	// Holds the information about the SUA / Filter
	FilteredAircraftDisplay fad = null;  // Not sure if this gets used anymore.
	
	List<FilteredAircraftDisplay> fadList = new ArrayList<FilteredAircraftDisplay>();
	List<String> fcaNameList = new ArrayList<String>();
	
	static boolean convertAsdiToTRX = false;
	static boolean filterRawData 	= false;
	static boolean simulateReroutes = false;
	static boolean simulateBaseline = true;

	
	public static void main(String[] args) {
		/*
		 * General stuff that'll be the same for every run
		 */
		String FACETdir = "/Users/tcolvin1/Documents/testWorkspace/testFacet";
		String baseDir = "/Users/tcolvin1/Documents/workspace3p6/facet/user/work/preferences/";
		String runName = "SimplePrefs";
		String prefskey = baseDir + runName;
		
//		for (String curArg : args){
//			System.out.println("Arg:" + curArg);
//		}
//		
//		System.exit(0);
		

		
		if (convertAsdiToTRX) {
			/*
			 * Run live data and dump to file.
			 * Look at the folders / datasets in the current folder.  Run them.
			 */		
			
			String dataHeadFolder = "/Users/tcolvin1/ToConvert/2013/";
			String outputFolder = "/Users/tcolvin1/Desktop/aiaaSpace/RawTRX/";
	
			// ==== Function ====
			//List<String> dataNameVector = new ArrayList();
			File yearFolder = new File(dataHeadFolder);		
				
			AiaaSpace curObj = new AiaaSpace();
			curObj.setup(FACETdir, prefskey);
	
			for (File monthFolder : yearFolder.listFiles()) {
				if (monthFolder.isDirectory()){   // Make sure that it's actually a folder
					//File[] listOfDays = monthFolder.listFiles();
					for (File dayFolder : monthFolder.listFiles()) {
						if (dayFolder.isDirectory()) {
							String outFileName = yearFolder.getName() + monthFolder.getName() + dayFolder.getName() + ".trx";
							System.out.println(outFileName);
											
							curObj.runLiveData(dayFolder.getAbsolutePath(), outputFolder + outFileName);
	
						}
					}
				}
			}
			
			System.exit(0);
		}
		

		
		
		if (filterRawData) {
			/*
			 * Run multiple filters
			 * BUG: If TRX files are out of date order, then this won't work, facet will simply advance the date one day.
			 *        So the traffic data will get labeled incorrectly.
			 */
			// String filterFolderName = "/Users/tcolvin1/Documents/SU-FARM/papers/Space2015/SpaceportFilters/";
			// String[] filterName = {"SUA_FrntRngeFilter", "SUA_MidlandFilter", "SUA_MojaveFilter"};
			
			String filterFolderName = "/Users/tcolvin1/Documents/SU-FARM/papers/Space2015/Filters/";		
			// String filterName = "SUA_Total";
			String filterName = "SUA_OKFilter";
			
			String dataFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/RawTRX/";
			File dataFolder = new File(dataFolderName);
			String outFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/FilteredTRX/";
			
	
			List<String> dataNameList = new ArrayList<String>();		
			String[] allFileList = dataFolder.list();
			for (String curFile : allFileList){
				if (curFile.endsWith(".trx")){
					dataNameList.add(curFile);
					System.out.println(curFile);
				}
			}
			
			System.out.println("RUNNING AiaaSpace from API");
			AiaaSpace curObj = new AiaaSpace();
			curObj.setup(FACETdir, prefskey);
			
			curObj.loadFilterBatch(filterFolderName, filterName);
			// curObj.loadSUA(filterFolderName + filterName, outFolderName + "test.flt");
	
			
			for (String dataName : dataNameList) {
				String data = dataFolderName + dataName;
				
				String baseDataName = dataName.substring(4, dataName.lastIndexOf(".trx"));
				System.out.println(data);
				
				curObj.updateOutfileNameBatch(outFolderName, "_" + baseDataName + ".flt");
				curObj.runFilter(data);		
			}
			
			System.exit(1);
		}
		
		
		
		if (simulateReroutes) {
			/*
			 * Trying to fix the memory leak.
			 * Run the space ops on all the downloaded days
			 * BUG: If TRX files are out of date order, then this won't work, facet will simply advance the date one day.
			 *        So the traffic data will get labeled incorrectly.
			 *        
			 */
	
			String outFolderName = args[0];
	//		String suaPathPlusFileName = args[1];
			File missionSuaFile = new File(args[1]);
	//		String dataName = args[2];
			File dataFile = new File(args[2]);
			
			String missionSuaFolderName = missionSuaFile.getParent() + "/";
			String curSUA = missionSuaFile.getName();
			
			String dataName = dataFile.getName();
			String dataFolderName = dataFile.getParent() + "/";
			
			
	//		System.out.println("suaFolder = " + missionSuaFolderName);
	//		System.out.println("curSUA    = " + curSUA);
			
			System.out.println(curSUA + "  " + dataName);
			
			String[] missionParts = curSUA.split("_");
			String scenario = missionParts[1];
			String vehicleName = missionParts[2];
			String launchLocation = missionParts[3];
			
			String baseDataName = dataName.substring(4, dataName.lastIndexOf(".trx"));
			String comment = scenario;
	
	//		System.exit(0);
			
	
			// Fire up the beast
			System.out.println("RUNNING AiaaSpace from API");
			AiaaSpace curObj = new AiaaSpace();
			curObj.setup(FACETdir, prefskey);
			
			// Load the SUA file (only one can be passed in)
			System.out.println("Loading: " + missionSuaFolderName + curSUA);
			curObj.loadSUA(missionSuaFolderName + curSUA);  // Loads the SUA file but delays naming the output file.
			
			// Only one data file allowed right now.  If more come in, then this becomes a loop.
			curObj.updateCustomFilterName(outFolderName,  vehicleName + "_" + launchLocation + "_" + baseDataName + "_" + comment + ".flt");
	
			curObj.createCustomOutfileName(outFolderName, vehicleName + "_" + launchLocation + "_" + baseDataName + "_" + comment + ".cst");
	
			//curObj.runFilter(data);	
			String data = dataFolderName + dataName;
			int lookaheadMinutes = 10;
			curObj.runAvoidanceBatch(data, lookaheadMinutes);
			
			// This is TOTALLY a hack, but it works.  If you accidentally call the create function twice 
			//  in a row, you'll be sorry.  This remove function simply deletes the last-created plot.  If
			//  that plot just happens to be the custom one you wanna axe, then all is well.
			curObj.removeCustomOutfileName(outFolderName, "_" + baseDataName + comment +".cst");  // These arguments are dummies.
			
			System.exit(0);
		}
		
		
		if (simulateBaseline) {
			/*
			 * Kevin - 2025 Baseline
			 */
			
			String outFolderName = args[0];
			File dataFile = new File(args[1]);
					
			String dataName = dataFile.getName();
			String baseDataName = dataName.substring(4, dataName.lastIndexOf(".trx"));
	
			String dataFolderName = dataFile.getParent() + "/";
			String data = dataFolderName + dataName;
	
			String suaFile = "";
	
			int lookaheadMinutes = 10;
			String customOutFile = outFolderName + baseDataName + ".cst";
			String outfilterName = outFolderName + baseDataName + ".flt";
			
	//		System.out.println(customOutFile);
	//		System.out.println("   " + outfilterName);
	//		System.exit(0);
			
			System.out.println("AIAASPACE: Running FACET from the API!");
			new AiaaSpace().doRun(FACETdir, prefskey, data, suaFile, customOutFile, outfilterName, lookaheadMinutes);
		}
		
	
		
		
		
		
		
		
		
		
		
//		
//		
//		// Where the SUA files to simulate are located
//		// Where to put the output .flt and .cst files
////		String missionSuaFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/SuaFilesBrokenOut/TraditionalSUAs/";
////		String outFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/TraditionalOneMonth/";
//
////		String missionSuaFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/SuaFilesBrokenOut/FacetEnvelopes/";
////		String outFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/EnvelopesOneMonth/";
//		
//		
//		// Where the post-filter post-prune TRX files are located
////		String dataFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/PrunedTRX/runTheseFirstMonth/";
//		
//		// Parse out the SUA_ files to simulate
//		File missionSuaFolder = new File(missionSuaFolderName);
//		String[] allFileList = missionSuaFolder.list();
//		
//		List<String> missionList = new ArrayList<String>();		
//		for (String curFile : allFileList){
//			if (curFile.startsWith("SUA_")){
//				missionList.add(curFile);
//				//System.out.println(curFile);
//			}
//		}
//		
//		// Parse out all of the possible *.trx files
//		File dataFolder = new File(dataFolderName);
//		allFileList = dataFolder.list();
//		
//		List<String> dataList = new ArrayList<String>();		
//		for (String curFile : allFileList){
//			if (curFile.endsWith(".trx")){
//				dataList.add(curFile);
//				//System.out.println(curFile);
//			}
//		}
//
////		// This returns the runtime object associated with the current java application
////		Runtime r = Runtime.getRuntime();
//
//		
////		// Fire up the beast
////		System.out.println("RUNNING AiaaSpace from API");
////		AiaaSpace curObj = new AiaaSpace();
////		curObj.setup(FACETdir, prefskey);
//				
//		// For each SUA_ file, find the filtered .trx files that correspond to the spaceport / launch location
//		for (String curSUA : missionList){
//			String[] missionParts = curSUA.split("_");
//			String scenario = missionParts[1];
//			String vehicleName = missionParts[2];
//			String launchLocation = missionParts[3];
//			
//			
//			// Skip all SUA files that aren't this one.
//			if (!curSUA.equals("SUA_2018H_Lynx_Midland")) {
//				continue;
//			}
//			
//			
//			// Find the data files for the current SUA_
//			List<String> matchedData = new ArrayList<String>();	
//			for (String curData : dataList) {
//				// Check if launch location is in data file name.  Skip the first 4 characters (TRX_)
//				if (curData.startsWith(launchLocation,4)) {
//					matchedData.add(curData);
//				}
//			}
//			
//			// Make sure you didn't goof on the spelling of a launch location or file name format
//			if (matchedData.size() == 0) {
//				System.err.println("ERROR: NO MATCHING DATA FILES WERE FOUND FOR " + curSUA);
//				System.exit(-1);
//			}
//
//			// Tell me what they are
//			System.out.println(curSUA + " will run:");
//			for (String curData : matchedData) {
//				System.out.println("   " + curData);				
//			}
//
////			// Run them!
////			System.out.println("Loading: " + missionSuaFolderName + curSUA);
////			curObj.loadSUA(missionSuaFolderName + curSUA);  // Loads the SUA file but delays naming the output file.
//
//			
//			for (String dataName : matchedData) {
//				// Fire up the beast
//				System.out.println("RUNNING AiaaSpace from API");
//				AiaaSpace curObj = new AiaaSpace();
//				curObj.setup(FACETdir, prefskey);
//				// Run them!
//				System.out.println("Loading: " + missionSuaFolderName + curSUA);
//				curObj.loadSUA(missionSuaFolderName + curSUA);  // Loads the SUA file but delays naming the output file.
//				
//				
//				String data = dataFolderName + dataName;
//				
//				String baseDataName = dataName.substring(4, dataName.lastIndexOf(".trx"));
//				// baseDataName looks like FrontRangeFilter_20130101
//				String[] parts = baseDataName.split("_");
//				String trafficDate = parts[1];
//				
//				System.out.println(data);
//				
//				String comment = scenario;
//
////				// Skip any input files before this date.
////				if (Integer.parseInt(trafficDate) < 20130123){
////					continue;
////				}
//				
//				curObj.updateCustomFilterName(outFolderName,  vehicleName + "_" + launchLocation + "_" + baseDataName + "_" + comment + ".flt");
//		
//				curObj.createCustomOutfileName(outFolderName, vehicleName + "_" + launchLocation + "_" + baseDataName + "_" + comment + ".cst");
//
//				//curObj.runFilter(data);	
//				int lookaheadMinutes = 10;
//				curObj.runAvoidanceBatch(data, lookaheadMinutes);
//				
//				// This is TOTALLY a hack, but it works.  If you accidentally call the create function twice 
//				//  in a row, you'll be sorry.  This remove function simply deletes the last-created plot.  If
//				//  that plot just happens to be the custom one you wanna axe, then all is well.
//				curObj.removeCustomOutfileName(outFolderName, "_" + baseDataName + comment +".cst");  // These arguments are dummies.
//				
//				curObj.tryStuff();
////				curObj.closeFacet();
//				
////				long availMem1 = r.freeMemory();
////				r.gc();
////				long availMem2 = r.freeMemory();
////				System.out.println("Suggested Garbage Collection Freed: " + (availMem1-availMem2) );
//				
//			}	
//		}
//		
//		System.exit(1);
//		
//		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	
		
//		/*
//		 * THIS WORKS!  But it has a memory leak or something so you can't run for too long.
//		 * Run the space ops on all the downloaded days
//		 * BUG: If TRX files are out of date order, then this won't work, facet will simply advance the date one day.
//		 *        So the traffic data will get labeled incorrectly.
//		 *        
//		 */
//		
////		String custom = "Custom_";
////		String filter = "Filter_";
////		
////		String cst = ".cst";
////		String flt = ".flt";
//
//		// Where the SUA files to simulate are located
//		// Where to put the output .flt and .cst files
////		String missionSuaFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/SuaFilesBrokenOut/TraditionalSUAs/";
////		String outFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/TraditionalOneMonth/";
//
//		String missionSuaFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/SuaFilesBrokenOut/FacetEnvelopes/";
//		String outFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/MissionResults/EnvelopesOneMonth/";
//
//		// Where the post-filter post-prune TRX files are located
//		String dataFolderName = "/Users/tcolvin1/Desktop/aiaaSpace/PrunedTRX/runTheseFirstMonth/";
//		
//		// Parse out the SUA_ files to simulate
//		File missionSuaFolder = new File(missionSuaFolderName);
//		String[] allFileList = missionSuaFolder.list();
//		
//		List<String> missionList = new ArrayList<String>();		
//		for (String curFile : allFileList){
//			if (curFile.startsWith("SUA_")){
//				missionList.add(curFile);
//				//System.out.println(curFile);
//			}
//		}
//		
//		// Parse out all of the possible *.trx files
//		File dataFolder = new File(dataFolderName);
//		allFileList = dataFolder.list();
//		
//		List<String> dataList = new ArrayList<String>();		
//		for (String curFile : allFileList){
//			if (curFile.endsWith(".trx")){
//				dataList.add(curFile);
//				//System.out.println(curFile);
//			}
//		}
//		
//		// Fire up the beast
//		System.out.println("RUNNING AiaaSpace from API");
//		AiaaSpace curObj = new AiaaSpace();
//		curObj.setup(FACETdir, prefskey);
//				
//		// For each SUA_ file, find the filtered .trx files that correspond to the spaceport / launch location
//		for (String curSUA : missionList){
//			String[] missionParts = curSUA.split("_");
//			String scenario = missionParts[1];
//			String vehicleName = missionParts[2];
//			String launchLocation = missionParts[3];
//			
//			
//			// Skip all SUA files that aren't this one.
//			if (!curSUA.equals("SUA_2018H_Lynx_Midland")) {
//				continue;
//			}
//			
//			
//			// Find the data files for the current SUA_
//			List<String> matchedData = new ArrayList<String>();	
//			for (String curData : dataList) {
//				// Check if launch location is in data file name.  Skip the first 4 characters (TRX_)
//				if (curData.startsWith(launchLocation,4)) {
//					matchedData.add(curData);
//				}
//			}
//			
//			// Make sure you didn't goof on the spelling of a launch location or file name format
//			if (matchedData.size() == 0) {
//				System.err.println("ERROR: NO MATCHING DATA FILES WERE FOUND FOR " + curSUA);
//				System.exit(-1);
//			}
//
//			// Tell me what they are
//			System.out.println(curSUA + " will run:");
//			for (String curData : matchedData) {
//				System.out.println("   " + curData);				
//			}
//
//			// Run them!
//			System.out.println("Loading: " + missionSuaFolderName + curSUA);
//			curObj.loadSUA(missionSuaFolderName + curSUA);  // Loads the SUA file but delays naming the output file.
//			//System.exit(-1);
//
//			for (String dataName : matchedData) {
//				String data = dataFolderName + dataName;
//				
//				String baseDataName = dataName.substring(4, dataName.lastIndexOf(".trx"));
//				// baseDataName looks like FrontRangeFilter_20130101
//				String[] parts = baseDataName.split("_");
//				String trafficDate = parts[1];
//				
//				System.out.println(data);
//				
//				String comment = scenario;
//
////				// Skip any input files before this date.
////				if (Integer.parseInt(trafficDate) < 20130123){
////					continue;
////				}
//				
//				curObj.updateCustomFilterName(outFolderName,  vehicleName + "_" + launchLocation + "_" + baseDataName + "_" + comment + ".flt");
//		
//				curObj.createCustomOutfileName(outFolderName, vehicleName + "_" + launchLocation + "_" + baseDataName + "_" + comment + ".cst");
//
//				//curObj.runFilter(data);	
//				int lookaheadMinutes = 10;
//				curObj.runAvoidanceBatch(data, lookaheadMinutes);
//				
//				// This is TOTALLY a hack, but it works.  If you accidentally call the create function twice 
//				//  in a row, you'll be sorry.  This remove function simply deletes the last-created plot.  If
//				//  that plot just happens to be the custom one you wanna axe, then all is well.
//				curObj.removeCustomOutfileName(outFolderName, "_" + baseDataName + comment +".cst");  // These arguments are dummies.
//			}	
//		}
//		
//		System.exit(1);
	

		
	
		

	}
	
	private void tryStuff() {
//		srv.getFACETRunInterface().restartSim()
//	    srv.getSim().terminate();
//	    srv.getFACETRunInterface().getFACETInterface().rmAllFlightHistories();
//	    srv.getFACETRunInterface().getFACETInterface().rmAllFlightPlans();
		
		srv.disposeFACET();
		srv = null;
		
		// This returns the runtime object associated with the current java application
		Runtime r = Runtime.getRuntime();
		r.gc();
		
	    

//		srv.disposeFACET();		
//		srv.getFACETRunInterface().getFACETInterface()
	}
	

	
	private void setup(String dir, String prefskey) {
		boolean runHeadless = false;
		
		//Initiates connection to FACET
		srv = null;
		try {
			srv = facet_pkg.api.server.FACETServerAPI.getInstance();	
			
			if (srv == null) {
				System.out.println("Couldn't make connection to FACET");
			} } 
		catch ( Exception e ) {e.printStackTrace(); System.exit(1); }
		
		// Opens FACET from given directory with preferences file
		File directory = new File(dir);
		File preferenceFile = new File(prefskey);
		try {
			srv.loadFACET(directory, preferenceFile, runHeadless);   // This should work!  Why is it commented out?
//			srv.loadFACET(directory, preferenceFile); 
			} 
		catch ( Exception e ) {
			System.out.println("Something went wrong with the loading preferences");
			e.printStackTrace(); 
			System.exit(1); }
		
//		srv.simStop()
//		srv.disposeFACET()
	}
	
	private void closeFacet() {
		srv.disposeFACET();		
//		srv.getFACETRunInterface().getFACETInterface()
	}

	private void loadSUA(String suaFile) {
		try {
		    /*
		     * Load the SUAs
		     */
			SUAInterface sua = srv.getSUA();
			if (suaFile.length() > 0) {
		        sua.loadUserSUAFile(suaFile);
			}
	
			int userType = sua.USER;
			int numSuas = sua.getCount(userType);
			String[] suaNameVec = new String[numSuas];
	        String fca = "";
	
			for (int ix = 0; ix < numSuas; ix++) {
				suaNameVec[ix] = sua.getName(userType, ix);
				fca = fca + suaNameVec[ix] + ",";
			}
			System.out.println("fca = " + fca);
			
			System.out.println("numSuas = " + numSuas);
			for (String curName : suaNameVec) {
				System.out.format("curName = %s \n", curName);
			}
	
			String acid = null;
			String regNum = null;
			String type = null;
	        String origin = null;
	        String destination = null;
	        String flight_level = null;
	        String center = null;
	        String sector = null;
	        String jet_route = null;
	        String origin_center = null;
	        String destination_center = null;
	        String fix = null;
	        String weight = null;
	        String equippage = null;
	        String haul = null;
	        String flight_plan = null;
	        String delay = null;
	//        String fca = ",".;
	        Boolean controlled = null;
	        String reroute = null;
	        String streamName = null;
	        
			AircraftFilterData filterData = new AircraftFilterData( acid, regNum, type,
		            origin, destination, flight_level,
		            center, sector, jet_route,
		            origin_center, destination_center,
		            fix, weight, equippage,
		            haul, flight_plan, delay,
		            fca,
		            controlled, reroute, streamName );
			
			AircraftFilterManager afm = srv.getFACETRunInterface().getEngineInterface().getAircraftFilterManager();
	
			// Clear out all previous filters and add the new one to the manager
			afm.removeAllAircraftFilters();
			afm.addAircraftFilter(filterData);
			
			
			// Then load the two display windows with the filter info
			AircraftFilterDialog afd = new AircraftFilterDialog(srv.getFACETRunInterface().getFACETInterface(), afm );
			fad = new FilteredAircraftDisplay(afm, filterData);
//			FilteredAircraftDisplay fad = new FilteredAircraftDisplay(afm, filterData);
			
			
			// Then connect the display we created with the data that we started with
			filterData.setFilteredAircraftDisplay(fad);
			
			// Tell the filter to dump continuously to a file.
			//		fad.useCustomFilterOutput = true;	//Tell it to use my custom file format (different from custom TRX format, haha)
//			fad.continuouslyWriteFile(outfilterName);
			
			// This is not a batch file, so only want a single fad.  Always clear the list before adding the current ONE.
			fadList.clear();
			fcaNameList.clear();
			fadList.add(fad);
			fcaNameList.add(fca);

		} finally {
//			System.out.println("SUAs were either not present or failed to load");
		}
	
	}

//	// This was a poorly-named function.  Change it to updateFilterName or something like that.
//	private void updateOutfileName(String outfilterName) {
//		fad.continuouslyWriteFile(outfilterName);
//	}

	private void updateCustomFilterName(String outFolderName, String comment) {
//		fad.continuouslyWriteFile(outfilterName);
		// Want to leave it general if possible, so use fadList
		// make sure fad.useCustomFilterOutput = true; is set somewhere too
		
		int curSize = fadList.size();
		if (curSize  > 1) {
			System.err.println("ERROR: updateCustomFilterName is not meant to be used with multiple fads in fadList!");
			System.exit(-1);
		} else if (curSize <= 0) {
			System.err.println("ERROR: updateCustomFilterName cannot accept a length-zero fadlist.  Must have exactly one element!");
			System.exit(-1);			
		}

		String outName = outFolderName + comment;
		fadList.get(0).continuouslyWriteFile(outName);
		
		
//		for (int ix = 0; ix < fcaNameList.size(); ix++){
//			//String outName = outFolderName + fcaNameList.get(ix) + comment;
//			String outName = outFolderName + comment;
//			fadList.get(ix).continuouslyWriteFile(outName);
//		}
		
	}
	

	
	private void createCustomOutfileName(String outFolderName, String customOutFile) {
//		fad.continuouslyWriteFile(outfilterName);
		
//		ERROR;
		// I don't think there's any (easy) way to get this to spit out multiple output files.
		// All I want is to just be able to update the one.
		
		/*
		 * This is a custom output file / format that I wrote which will give us everything except for the 
		 * 		flight plans of each flight.  If any changes are desired, edit writeCustomOutputTJC.c
		 */
		
		// Heck, I may not even have to mess with the fadList?
		GUIInterface gi = srv.getGUI();
		gi.recordCustomTJC(outFolderName + customOutFile);
	}
	
	private void removeCustomOutfileName(String outFolderName, String customOutFile) {
		GUIInterface gi = srv.getGUI();
		gi.removeCustomTJC("balls");
//		removeCustomTJC
	}
	
	private void loadFilterBatch(String filterFolderName, String filterName) {
		/*
		 * Pass in a filter (SUA_ file). This will take each individual shape defined within
		 * 		and associate with it its own output file.  
		 * 		+ Great for lots of independent shapes.
		 * 		+ Bad for dynamic / related shapes.  
		 * 		+ See loadSUABatch for similar but different functionality.
		 */
		
		try {
		    /*
		     * Load the SUAs
		     */

			SUAInterface sua = srv.getSUA();
	
			String suaFile = filterFolderName + filterName;
			if (suaFile.length() > 0) {
		        sua.loadUserSUAFile(suaFile);
			}
					
			int userType = sua.USER;
			int numSuas = sua.getCount(userType);
			String[] suaNameVec = new String[numSuas];
//	        String fca = "";
//	
			for (int ix = 0; ix < numSuas; ix++) {
				suaNameVec[ix] = sua.getName(userType, ix);
//				fca = fca + suaNameVec[ix] + ",";
			}
//			System.out.println("fca = " + fca);
//			
//			System.out.println("numSuas = " + numSuas);
//			for (String curName : suaNameVec) {
//				System.out.format("curName = %s \n", curName);
//			}
			
			for (int ix = 0; ix < numSuas; ix++) {
		        String fca = suaNameVec[ix];  // Each SUA will be a single filter
		        System.out.println("fca = " + fca);

				
				String acid = null;
				String regNum = null;
				String type = null;
		        String origin = null;
		        String destination = null;
		        String flight_level = null;
		        String center = null;
		        String sector = null;
		        String jet_route = null;
		        String origin_center = null;
		        String destination_center = null;
		        String fix = null;
		        String weight = null;
		        String equippage = null;
		        String haul = null;
		        String flight_plan = null;
		        String delay = null;
		//        String fca = ",".;
		        Boolean controlled = null;
		        String reroute = null;
		        String streamName = null;
		        
				AircraftFilterData filterData = new AircraftFilterData( acid, regNum, type,
			            origin, destination, flight_level,
			            center, sector, jet_route,
			            origin_center, destination_center,
			            fix, weight, equippage,
			            haul, flight_plan, delay,
			            fca,
			            controlled, reroute, streamName );
				
				AircraftFilterManager afm = srv.getFACETRunInterface().getEngineInterface().getAircraftFilterManager();
		
				// Add the filter to the manager
				afm.addAircraftFilter(filterData);
				
				// Then load the two display windows with the filter info
				AircraftFilterDialog afd = new AircraftFilterDialog(srv.getFACETRunInterface().getFACETInterface(), afm );
				FilteredAircraftDisplay curFad = new FilteredAircraftDisplay(afm, filterData);

//				fad = new FilteredAircraftDisplay(afm, filterData);
	//			FilteredAircraftDisplay fad = new FilteredAircraftDisplay(afm, filterData);
				
				// Then connect the display we created with the data that we started with
				filterData.setFilteredAircraftDisplay(curFad);
				
				// Tell the filter to dump continuously to a file.
				//		fad.useCustomFilterOutput = true;	//Tell it to use my custom file format (different from custom TRX format, haha)
//				fad.continuouslyWriteFile(outfilterName);
				
				fadList.add(curFad);
				fcaNameList.add(fca);
			}
			
		} finally {
//			System.out.println("SUAs were either not present or failed to load");
		}
	
	}

	
	private void loadSUABatch(String filterFolderName, String filterName) {
		
		try {
		    /*
		     * Load the SUAs
		     */

			SUAInterface sua = srv.getSUA();
	
			String suaFile = filterFolderName + filterName;
			if (suaFile.length() > 0) {
		        sua.loadUserSUAFile(suaFile);
			}
					
			int userType = sua.USER;
			int numSuas = sua.getCount(userType);
			String[] suaNameVec = new String[numSuas];
//	        String fca = "";
//	
			for (int ix = 0; ix < numSuas; ix++) {
				suaNameVec[ix] = sua.getName(userType, ix);
//				fca = fca + suaNameVec[ix] + ",";
			}
//			System.out.println("fca = " + fca);
//			
//			System.out.println("numSuas = " + numSuas);
//			for (String curName : suaNameVec) {
//				System.out.format("curName = %s \n", curName);
//			}
			
			for (int ix = 0; ix < numSuas; ix++) {
		        String fca = suaNameVec[ix];  // Each SUA will be a single filter
		        System.out.println("fca = " + fca);

				
				String acid = null;
				String regNum = null;
				String type = null;
		        String origin = null;
		        String destination = null;
		        String flight_level = null;
		        String center = null;
		        String sector = null;
		        String jet_route = null;
		        String origin_center = null;
		        String destination_center = null;
		        String fix = null;
		        String weight = null;
		        String equippage = null;
		        String haul = null;
		        String flight_plan = null;
		        String delay = null;
		//        String fca = ",".;
		        Boolean controlled = null;
		        String reroute = null;
		        String streamName = null;
		        
				AircraftFilterData filterData = new AircraftFilterData( acid, regNum, type,
			            origin, destination, flight_level,
			            center, sector, jet_route,
			            origin_center, destination_center,
			            fix, weight, equippage,
			            haul, flight_plan, delay,
			            fca,
			            controlled, reroute, streamName );
				
				AircraftFilterManager afm = srv.getFACETRunInterface().getEngineInterface().getAircraftFilterManager();
		
				// Add the filter to the manager
				afm.addAircraftFilter(filterData);
				
				// Then load the two display windows with the filter info
				AircraftFilterDialog afd = new AircraftFilterDialog(srv.getFACETRunInterface().getFACETInterface(), afm );
				FilteredAircraftDisplay curFad = new FilteredAircraftDisplay(afm, filterData);

//				fad = new FilteredAircraftDisplay(afm, filterData);
	//			FilteredAircraftDisplay fad = new FilteredAircraftDisplay(afm, filterData);

				// Not 100% sure, but I think this needs to be set BEFORE we add it to the filtered AC display.
				curFad.useCustomFilterOutput = true;  // Hopefully this works!

				// Then connect the display we created with the data that we started with
				filterData.setFilteredAircraftDisplay(curFad);
				
				// Tell the filter to dump continuously to a file.
				//		fad.useCustomFilterOutput = true;	//Tell it to use my custom file format (different from custom TRX format, haha)
//				fad.continuouslyWriteFile(outfilterName);
				
				
				fadList.add(curFad);
				fcaNameList.add(fca);
			}
			
		} finally {
//			System.out.println("SUAs were either not present or failed to load");
		}
	
	}
	
	private void updateOutfileNameBatch(String outfilterFolderName, String comment) {
		if (fcaNameList.size() != fadList.size()){
			System.out.println("updateOutfileName failed");
			System.exit(-2);
		}
		
		for (int ix = 0; ix < fcaNameList.size(); ix++){
			String outName = outfilterFolderName + fcaNameList.get(ix) + comment;
			fadList.get(ix).continuouslyWriteFile(outName);
		}
		
//		fadList.set(ix, element);
//		fadList.get(ix).continuouslyWriteFile(outfilterName);
	}
	
	

	
	private void runFilter(String inputTRXFile) {
        /*
         * Load up the simulation
         */
	    int trajectoryUpdateInterval 	= 60;	// [sec] How often data gets output / screen displays.
	    double integrationTimeStep 		= 5.0;	// [sec] Presumably this is self-explanatory.
	    double additionalUpdateDelay 	= 0.0;	// [sec] Used to SLOW DOWN the GUI for visualization purposes.
        	    
		srv.getSim().startSimOnPause(true);
		srv.getSim().loadPlaybackAsynch(inputTRXFile, trajectoryUpdateInterval, additionalUpdateDelay);
		
//		/*
//		 * Run the simulation and close FACET when it finishes
//		 */
		srv.getSim().resume();
//		srv.getSim().pauseResume();
		
		/*
		 * Kill FACET once the simulation completes
		 */
		try {
			synchronized (srv) {
				while (srv.getSim().isRunning()) {
				srv.wait(2000);	// Wait 2 seconds
				} }
			System.out.println("isRunning = " + srv.getSim().isRunning());
//			srv.getSim().pauseResume();
			
			srv.getSim().terminate();


//			System.exit(0);
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
	}
	
	
	
	private void runLiveData(String input_directory_str, String trxOutFile) {
        /*
         * Load up the simulation
         */
	    int trajectoryUpdateInterval 	= 60;	// [sec] How often data gets output / screen displays.
	    double integrationTimeStep 		= 5.0;	// [sec] Presumably this is self-explanatory.
	    double additionalUpdateDelay 	= 0.0;	// [sec] Used to SLOW DOWN the GUI for visualization purposes.
        	    
//		srv.getSim().startSimOnPause(true);
//		srv.getSim().loadPlaybackAsynch(inputTRXFile, trajectoryUpdateInterval, additionalUpdateDelay);
		
		
//		/*
//		 * Run the simulation and close FACET when it finishes
//		 */
//		srv.getSim().resume();
//		srv.getSim().pauseResume();
	    
	    
//		/*
//		 * Save a TRX file of what was simulated.  This will have all flight plans / changes inside
//		 */
//	    long secondsSinceMidnightStart 	= startTime % (24*3600);
//	    long secondsSinceMidnightEnd 	= (startend[1]/1000) - secondsSinceMidnightStart + 24*3600;	//Record for one full day past end of TRX
////	    long secondsSinceMidnightEnd 	= (startend[1]/1000) % (24*3600);
//	    
		// String outputTRX = "TRX_LevelOneOutput";
	    
	    
	    
	    long secondsSinceMidnightStart = 1;	            // 0 silently fails for some reason
	    long secondsSinceMidnightEnd = 365*24*60*60;	//Not sure how long files are, so basically don't ever stop
	    
	    RcdAcPosSetup curSetup = new RcdAcPosSetup(trxOutFile, secondsSinceMidnightStart, secondsSinceMidnightEnd, true);
		srv.getFACETRunInterface().getFACETInterface().setRcdAcPosSetup(curSetup);  
	    
//		srv.getFACETRunInterface().getFACETInterface().
//			setRcdAcPosSetup(new RcdAcPosSetup(trxOutFile, secondsSinceMidnightStart, secondsSinceMidnightEnd, true));   
	   
		    
		/*
		 * Look at the files in the current folder.  Figure out which ones will trigger the facet timestep error
		 * 		and just delete them. Have to get rid of the first 5 minutes OR SO.  Some days need the first 6...
		 */	

		File folder = new File(input_directory_str);
		File[] listOfFiles = folder.listFiles();
		int numNonDataFiles = 0;

		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile()) {
				String curName = listOfFiles[i].getName();
				if (curName.contains(".xml.bin")){
					int timeStr = Integer.parseInt( curName.substring(17, 21) );
					if (timeStr <= 5) {  
						System.out.println("DELETING file " + curName);
						listOfFiles[i].delete();
					}
				} else {
					// Keep track of how many other files are in the folder so that we know when we're done.
					numNonDataFiles++;
				}
			}
		}
		
		    
	   
	    
	    // The ac_list file is weird and I dunno what it is.  It's definitely not what I want.  Don't let it happen.
	    //   I think it also has a bug, because the file size was *decreasing* for awhile...so it must not append.
	    //   Unfortunately, it seems like you have to provide a valid file / folder, otherwise FACET will hang.
//		String testStr = "/Users/tcolvin1/Desktop/Jan 2013/DW2015013118570152/";
		String preprocessor_input_file_str = null;
//		String save_ac_list_file_str = "/Users/tcolvin1/Desktop/aiaaSpace/Results";
		String save_ac_list_file_str = "trashFile";
		String save_ac_list_freq_str = "0";  // This frequency is timesteps, not seconds.  I think zero means don't ever write.
		String load_ac_list_file_str = null;
		
		String filename_prefix = "asdi_xml.";
		
		System.out.println("input_directory_str: " + input_directory_str);
		
		srv.getSim().runLiveASDIXML(true, 0.05, filename_prefix, input_directory_str, 
				preprocessor_input_file_str, save_ac_list_file_str, save_ac_list_freq_str, load_ac_list_file_str);
		

		
		/*
		 * Kill FACET once the simulation completes
		 */
		try {
			synchronized (srv) {
				// Run until all of the data files have been used up
				while (folder.listFiles().length != numNonDataFiles ) {
//					while (srv.getSim().isRunning()) {
//					listOfFiles = folder.listFiles();
//					listOfFiles.length != numNonDataFiles;
					srv.wait(2000);	// Wait 2 seconds
				} 
				srv.wait(10000);	// Wait 10 more seconds just to make sure everything is run through
				
			}
//			srv.getSim().pauseResume();
			srv.getFACETRunInterface().getFACETInterface().stopCurrentRun();
			srv.getSim().terminate();

			System.out.println("isRunning = " + srv.getSim().isRunning());

//			System.exit(0);
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
	}

	
	private void runAvoidanceBatch(String inputTRXFile, int lookaheadMinutes) {
        /*
         * Load up the simulation
         */
	    int trajectoryUpdateInterval 	= 60;	// [sec] How often data gets output / screen displays.
	    double integrationTimeStep 		= 5.0;	// [sec] Presumably this is self-explanatory.
	    double additionalUpdateDelay 	= 0.0;	// [sec] Used to SLOW DOWN the GUI for visualization purposes.
        	    
		srv.getSim().startSimOnPause(true);
		srv.getSim().loadFlightPlanSimAsynch(inputTRXFile, trajectoryUpdateInterval, integrationTimeStep, additionalUpdateDelay);
		srv.getGUI().setAircraftSymbol(0);	// Triangle...actually you could take care of this with a preference file

		//		srv.getSim().loadPlaybackAsynch(inputTRXFile, trajectoryUpdateInterval, additionalUpdateDelay);
		
        Avoidance avoidance = srv.getFACETRunInterface().getFACETInterface().getAvoidance();
		if (fadList.size() > 0) {	        
	        // Start Avoidance Testing
			avoidance.avoidFCA(true);
			avoidance.lookAheadTime(lookaheadMinutes*60);  //10 min
			avoidance.setAlgorithm(AvoidanceAlgorithm.ConvexHull);
//			avoidance.setAlgorithm(AvoidanceAlgorithm.MinimumPath);
			avoidance.buffer(5.);
			avoidance.minimumFCASize(0.);
			avoidance.setupCellDetection();			
		}
		

//		/*
//		 * Run the simulation and close FACET when it finishes
//		 */
		srv.getSim().resume();
//		srv.getSim().pauseResume();
		
		/*
		 * Kill FACET once the simulation completes
		 */
		try {
			synchronized (srv) {
				while (srv.getSim().isRunning()) {
				srv.wait(2000);	// Wait 2 seconds
				} }
			System.out.println("isRunning = " + srv.getSim().isRunning());
//			srv.getSim().pauseResume();
			
			srv.getSim().terminate();


//			System.exit(0);
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
	}


	
	
	private void doRun(String dir, String prefskey, String inputTRXFile, String suaFile, String customOutFile, String outfilterName, int lookaheadMinutes) {
		boolean runHeadless = false;
		
		//Initiates connection to FACET
		FACETServerAPI srv = null;
		try {
			srv = facet_pkg.api.server.FACETServerAPI.getInstance();	
			
			if (srv == null) {
				System.out.println("Couldn't make connection to FACET");
			} } 
		catch ( Exception e ) {e.printStackTrace(); System.exit(1); }
		
		// Opens FACET from given directory with preferences file
		File directory = new File(dir);
		File preferenceFile = new File(prefskey);
		try {
			srv.loadFACET(directory, preferenceFile, runHeadless); } 
		catch ( Exception e ) {
			System.out.println("Something went wrong with the loading preferences");
			e.printStackTrace(); 
			System.exit(1); }		
		
		/* 
		 * The output will update every trajectoryUpdateInterval seconds (probably 60)
		 * 		and it would be nicest if these updates fell on a whole minute.  This
		 * 		checks the start time of the TRX file and makes sure that it falls nicely
		 * 		on a whole minute.
		 */
	    long[] startend = srv.getSim().getBeginEndTimeStamp(inputTRXFile);
	    
	    long startTime = startend[0]/1000;
	    if ((startTime % 60) != 0){
	    	System.err.println("ERROR: TRX file doesn't start on an even minute.  " +
	    			"Prepend TRACK_TIME " + (startTime - (startTime % 60)));
	    	System.exit(-1);
	    }
        
	    	    
        /*
         * Load up the simulation
         */
	    int trajectoryUpdateInterval 	= 60;	// [sec] How often data gets output / screen displays.
	    double integrationTimeStep 		= 5.0;	// [sec] Presumably this is self-explanatory.
	    double additionalUpdateDelay 	= 0.0;	// [sec] Used to SLOW DOWN the GUI for visualization purposes.
        	    
		srv.getSim().startSimOnPause(true);
		srv.getSim().loadFlightPlanSimAsynch(inputTRXFile, trajectoryUpdateInterval, integrationTimeStep, additionalUpdateDelay);
		srv.getGUI().setAircraftSymbol(0);	// Triangle...actually you could take care of this with a preference file
		
		
        SUAInterface sua = srv.getSUA();
        Avoidance avoidance = srv.getFACETRunInterface().getFACETInterface().getAvoidance();

		if (suaFile.length() > 0) {
		    /*
		     * Load the SUAs
		     */
	        sua.loadUserSUAFile(suaFile);
	        
	        // Start Avoidance Testing
			avoidance.avoidFCA(true);
			avoidance.lookAheadTime(lookaheadMinutes*60);  //10 min
			avoidance.setAlgorithm(AvoidanceAlgorithm.ConvexHull);
//			avoidance.setAlgorithm(AvoidanceAlgorithm.MinimumPath);
			avoidance.buffer(5.);
			avoidance.minimumFCASize(0.);
			avoidance.setupCellDetection();			
		}

		int userType = sua.USER;
		int numSuas = sua.getCount(userType);
		String[] suaNameVec = new String[numSuas];
        String fca = "";

		for (int ix = 0; ix < numSuas; ix++) {
			suaNameVec[ix] = sua.getName(userType, ix);
			fca = fca + suaNameVec[ix] + ",";
		}
		System.out.println("fca = " + fca);
		
		System.out.println("numSuas = " + numSuas);
		for (String curName : suaNameVec) {
			System.out.format("curName = %s \n", curName);
		}

		
		
		String acid = null;
		String regNum = null;
		String type = null;
        String origin = null;
        String destination = null;
        String flight_level = null;
        String center = null;
        String sector = null;
        String jet_route = null;
        String origin_center = null;
        String destination_center = null;
        String fix = null;
        String weight = null;
        String equippage = null;
        String haul = null;
        String flight_plan = null;
        String delay = null;
//        String fca = ",".;
        Boolean controlled = null;
        String reroute = null;
        String streamName = null;
        
		AircraftFilterData filterData = new AircraftFilterData( acid, regNum, type,
	            origin, destination, flight_level,
	            center, sector, jet_route,
	            origin_center, destination_center,
	            fix, weight, equippage,
	            haul, flight_plan, delay,
	            fca,
	            controlled, reroute, streamName );
		
		AircraftFilterManager afm = srv.getFACETRunInterface().getEngineInterface().getAircraftFilterManager();
//		String outfilterName = "/Users/tcolvin1/Desktop/KevinOutputs/testFilter";

		// Add the filter to the manager
		afm.addAircraftFilter(filterData);
		
		// Then load the two display windows with the filter info
		AircraftFilterDialog afd = new AircraftFilterDialog(srv.getFACETRunInterface().getFACETInterface(), afm );
		FilteredAircraftDisplay fad = new FilteredAircraftDisplay(afm, filterData);

		// Tell the filter to dump continuously to a file.
		fad.useCustomFilterOutput = true;	//Tell it to use my custom file format (different from custom TRX format, haha)
		fad.continuouslyWriteFile(outfilterName);
		
		// Then connect the display we created with the data that we started with
		filterData.setFilteredAircraftDisplay(fad);
		
		/*
		 * This is a custom output file / format that I wrote which will give us everything except for the 
		 * 		flight plans of each flight.  If any changes are desired, edit writeCustomOutputTJC.c
		 */
		GUIInterface gi = srv.getGUI();
		gi.recordCustomTJC(customOutFile);
		
		
		/*
		 * Run the simulation and close FACET when it finishes
		 */
		srv.getSim().resume();
		
		/*
		 * Kill FACET once the simulation completes
		 */
		try {
			synchronized (srv) {
				while (srv.getSim().isRunning()) {
				srv.wait(2000);	// Wait 2 seconds
				} }
			System.out.println("isRunning = " + srv.getSim().isRunning());
			System.exit(0);
		
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
	}	
	
	
}













