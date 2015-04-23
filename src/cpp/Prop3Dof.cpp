// Idea: Propogate everything as modular as possible.  So do all of the first stage and dump to file.  
//	Probably no need to do second stage.  Give every trajectory an ID number.  Don't worry about having
//  same weather profiles for different parts of the same trajectory (first stage up, first stage down, debris, etc).
//  Include timing uncertainty after propagations complete by just shifting start times around.

// I should have used BOOST! :(
// www.boost.org/doc/libs/1_50_0/

#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#include <map> // header file needed for to use MAP STL
using std::map;

#include "Trajectory.h"
#include "Debris.h"
#include "Architecture.h"
#include "Footprint3D.h"
#include "SkyGrid.h"
#include "PointCloud.h"

#include "timer.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#define Do_Propagate 1
#define Do_Stitch 2
#define Do_MakeTube 3
#define Do_FootprintFirstPass 4
#define Do_MakeFacetFiles 5
#define Do_MakeGoogleEarthFiles 6
#define Do_SimLeftRightHotCold 7
#define Do_MC 8
#define Merge_Footprint_Vectors 9

#define Do_GenerateSimpleIndividualDebrisFootprints 100
#define Do_CombineSimpleIndividualDebrisFootprints  101

#define Do_MakeGoogleEarthNominalTrajectory         200
#define Do_SANDBOX                                  300


int main (int argc, char * const argv[]) {
    
	int WhatToDo = Do_SANDBOX;
    
    // Set Atmosphere Options
//    char StateOption[] = "full rocket";
	char WindOption[] = "simple atmosphere";
	char DensityOption[] = "cantwell density";
//	string cape_nominal_traj_filename = "Files/zzz_nominal_cape_finished.txt";
	string cape_nominal_traj_filename = "Files/Falcon9CapePropOpt7km.txt";
	string mars_nominal_traj_filename = "Files/zzz_nominal_MARS.txt";
    
//	string nominal_traj_filename = "Files/Falcon9Cape15Clamped.txt";        //This actually kinda travels inland, debris reaches central florida
    double delta_t = 2.0;

    string nominalFileName("GeneratedFiles/trajectory_points.dat");
//	char OutputTrajectoryFile[] = "GeneratedFiles/mission1.dat";
    
	char CapeLrhcFile[] = "GeneratedFiles/CapeLRHC.dat";
	char MarsLrhcFile[] = "GeneratedFiles/MarsLRHC.dat";

	char CapeStageDownFile[] = "GeneratedFiles/CapeStageDown.dat";
	char MarsStageDownFile[] = "GeneratedFiles/MarsStageDown.dat";
    
    char CapeGoogleEarthFile[] = "GeneratedFiles/CapeFootPrint.kml";
    char MarsGoogleEarthFile[] = "GeneratedFiles/MarsFootPrint.kml";
    char FullMissionGoogleEarthFile[] = "GeneratedFiles/FullMissionFootPrint.kml";

    char CapeGoogleEarthTrajFile[] = "GeneratedFiles/CapeTraj.kml";
    char MarsGoogleEarthTrajFile[] = "GeneratedFiles/MarsTraj.kml";
    char Ss2GoogleEarthTrajFile[] = "GeneratedFiles/SS2Traj.kml";
    
    char TestSkyGridGoogleEarthFile[] = "GeneratedFiles/TestSkyGrid.kml";
    
    string CapeFootprintVectorFile("GeneratedFiles/CapeFootprintVector_Upto50km.dat");
    string MarsFootprintVectorFile("GeneratedFiles/MarsFootprintVector.dat");

    double tstepMinutes = 1.;

//    double binSizeKm = 5;
    int numPerBatch = 10;
    int numToWrite = numPerBatch;
    
    double zBinHeight = 5.;
    double xBinLength = 2;
    double yBinLength = 2;
    
    switch (WhatToDo) {
            
            
            
            
            
            
            
        case Do_SANDBOX:{
            
//            enum color_t {blue, black, white, red, purple };
//            const char *types[] = {'blue', 'black', 'white', 'red', 'purple'};
//            
//            color_t myColor;
//            
//            map< string, color_t>myMap;

            

        } break;
            
            
            
            
            
            
            
            
            
            
            
        case Do_MakeGoogleEarthNominalTrajectory:{
            char StateOption[] = "full rocket";
            Trajectory MyCape(StateOption, WindOption, DensityOption);
            MyCape.Initialize_Stages(numPerBatch, delta_t, cape_nominal_traj_filename);  //seed is still fixed!!!
            MyCape.PropagateWholeRocket();
            int timesteupsUntilCutoff = MyCape.getTimestepsUntilCutoff();
            MyCape.write_to_google_earth_native(CapeGoogleEarthTrajFile, numToWrite);   // Makes a google earth movie of TRAJECTORIES, not compact
            
            cout << "cutoff in seconds = " << timesteupsUntilCutoff*delta_t << endl;
            
        } break;
            
            
            
            
            
            
            
        case Do_GenerateSimpleIndividualDebrisFootprints:{
            // ~~~~~~~~~~~~~~~  Create a bunch of individual debris footprints ~~~~~~~~~~~~~~~~~~~~~~~
            // Generate a bunch of footprints at different explosion times
            string CapeGoogleEarthFileITER("GeneratedFiles/CapeFootPrint");
            string datExtension(".dat");
            string kmlExtension(".kml");
            
            
            for (double timeExplodeUntil = 90; timeExplodeUntil <= 90.; timeExplodeUntil += 5.){
                char datfname[CapeGoogleEarthFileITER.size() + 2 + 4];
                sprintf(datfname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, datExtension.c_str() );
                
                char kmlfname[CapeGoogleEarthFileITER.size() + 2 + 4];
                sprintf(kmlfname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, kmlExtension.c_str() );
                
                //        sprintf(fname, "%s%d%s",  altExplodeUntil, extension.c_str() );
                cout << "datfname = " << datfname << endl;
                cout << "kmlfname = " << kmlfname << endl;
                
                {
                    char StateOption[] = "full rocket";
                    Trajectory MyCape(StateOption, WindOption, DensityOption);
                    MyCape.Initialize_Stages(numPerBatch, delta_t, cape_nominal_traj_filename);  //seed is still fixed!!!
                    MyCape.PropagateWholeRocket();
                    int timesteupsUntilCutoff = MyCape.getTimestepsUntilCutoff();
                    MyCape.write_to_google_earth_native(CapeGoogleEarthTrajFile, numToWrite);   // Makes a google earth movie of TRAJECTORIES, not compact envelope
                    MyCape.write_LRHC_points_file(CapeLrhcFile, tstepMinutes);              // Exports rocket launch points as lat lon for making footprint
                    MyCape.write_stage_down_points_file(CapeStageDownFile, tstepMinutes);
                    
                    Footprint3D MyCapeFootprint(CapeLrhcFile, zBinHeight);
                    MyCapeFootprint.RemoveFootprintKeepTiming();
                    MyCape.MonteCarloDebris(MyCapeFootprint, timeExplodeUntil);
                    
                    MyCape.HowLongUntilEnterNAS();
                    cout << "storing as vector~~~~~~~~~~~~~~" << endl;
                    MyCapeFootprint.store_footprint_as_vector(datfname);
                    //                cout << "exporting~~~~~~~~~~~~~~~~~~~~" << endl;
                    MyCapeFootprint.exportGoogleEarth(kmlfname);
                    
                } }
            
            // ~~~~~~~~~~~~~~~~ End individual debris clouds genreation ~~~~~~~~~~~~~~~~~~~~~
        } break;

            
            
        case Do_CombineSimpleIndividualDebrisFootprints:{
            // ~~~~~~~~~~~ Create an updating footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            Footprint3D updatingFootprintZeroToOne("GeneratedFiles/CapeFootPrint 10.dat");
            
            // Generate a bunch of footprints at different explosion times
            string CapeGoogleEarthFileITER("GeneratedFiles/CapeFootPrint");
            string extension(".dat");
            
            for (double timeExplodeUntil = 20; timeExplodeUntil < 131.; timeExplodeUntil += 10.){
                char fname[CapeGoogleEarthFileITER.size() + 3 + 4];
                sprintf(fname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, extension.c_str() );
                cout << "fname = " << fname << endl;
                
                Footprint3D tempFootprint(fname);
                updatingFootprintZeroToOne.MergeFootprintVectors(tempFootprint);
                
            }
            cout << "launcLatTEMP = " << updatingFootprintZeroToOne.getLaunchLat() << endl;
            
            updatingFootprintZeroToOne.ProjectAllPointsDown();
            
            cout << "launcLatTEMP = " << updatingFootprintZeroToOne.getLaunchLat() << endl;
            
            updatingFootprintZeroToOne.ChopTimeAt(1);
            
            
            Footprint3D updatingFootprintOneToTwo("GeneratedFiles/CapeFootPrint 60.dat");
            for (double timeExplodeUntil = 70; timeExplodeUntil < 131.; timeExplodeUntil += 10.){
                char fname[CapeGoogleEarthFileITER.size() + 3 + 4];
                sprintf(fname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, extension.c_str() );
                cout << "fname = " << fname << endl;
                
                Footprint3D tempFootprint(fname);
                updatingFootprintOneToTwo.MergeFootprintVectors(tempFootprint);
                
            }
            updatingFootprintOneToTwo.ProjectAllPointsDown();
            updatingFootprintOneToTwo.ChopTimeAt(1);    //makes this fp last for one minute
            updatingFootprintOneToTwo.AddToFootprintUTC(0, 1);  //plus hours,minutes
            
            Footprint3D updatingFootprint;
            updatingFootprint = updatingFootprintZeroToOne;
            updatingFootprint.MergeFootprintVectors(updatingFootprintOneToTwo);
            
            // Translate footprint to Georgia  30.921463,-81.513451
            updatingFootprint.ChangeLaunchSiteToDeg(30.921463, -81.513451);
            updatingFootprint.exportGoogleEarth(CapeGoogleEarthFile);
            
            char CapeGoogleEarthFile2[] = "GeneratedFiles/CapeFootPrintBlowStaging.kml";
            Footprint3D FootprintStaging("GeneratedFiles/CapeFootPrint170.dat");
            FootprintStaging.ChangeLaunchSiteToDeg(30.921463, -81.513451);
            FootprintStaging.exportGoogleEarth(CapeGoogleEarthFile2);
            
            
            char CapeGoogleEarthFile3[] = "GeneratedFiles/CapeFootPrintBlowMaxQ.kml";
            Footprint3D FootprintMaxQ("GeneratedFiles/CapeFootPrint 90.dat");
            FootprintMaxQ.exportGoogleEarth(CapeGoogleEarthFile3);
            
            char CapeGoogleEarthFile4[] = "GeneratedFiles/CapeFootPrintLRHC.kml";
            Footprint3D FootprintLHRC(CapeLrhcFile, zBinHeight);
            Footprint3D FootprintEarlyDebris("GeneratedFiles/CapeFootPrint 70.dat");
            FootprintLHRC.MergeFootprintVectors(FootprintEarlyDebris);
            FootprintLHRC.exportGoogleEarth(CapeGoogleEarthFile4);
            
            // ~~~~~~~~~~~~~~ End create updating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            Footprint3D TotalMission;
            //        TotalMission = updatingFootprint;
            TotalMission = FootprintStaging;
            //        TotalMission = MarsFootprint;
            //        MyFootprint.store_footprint_as_points();
            //        TotalMission.exportGoogleEarth(FullMissionGoogleEarthFile);
            
            // Make FACET files
            bool makeFacet = true;
            if (makeFacet) {
                int launchTimeHours = 9;
                int launchTimeMinutes = 45;
                int offsetTimeMinutes = 0;	//This turns on the first SUA $offset minutes earlier than the launch time
                
                
                string folderPath("GeneratedFacet/");
                string folderBaseName("Sandbox");
                char buffer[40];
                sprintf(buffer,"__Time_%02d%02d", (int) floor(launchTimeHours),(int) floor(launchTimeMinutes));
                string timeStartInfo(buffer);
                
                sprintf(buffer, "__dTime_%02d", (int) floor(tstepMinutes));
                string timeDeltaInfo(buffer);
                
                sprintf(buffer,"__BinSizeKm_%.2f", zBinHeight);
                string binInfo(buffer);
                
                string folderName = folderPath + folderBaseName + timeStartInfo + timeDeltaInfo + binInfo + "/";
                
                string mkdir("mkdir ");
                system((mkdir + folderName).c_str());
                
                cout << (mkdir + folderName).c_str() << endl;
                
                TotalMission.make_facet_files(folderName, launchTimeHours*60 + launchTimeMinutes, offsetTimeMinutes, tstepMinutes);
            }

        } break;
            
            
            
            
            
        default:{
			cout << "ERROR!  You chose some WhatToDo option that doesn't exist ~~~~~~~~~~~~~~~~~~" << endl;
        } break;

            
            
            
    }
    
    string endMessage("say \"finished\"");
	system(endMessage.c_str());
    
    return 0;
}
    

    
    
    
    
//    bool fresh = false;
    
    
    
//    if (fresh){
    

//        // ~~~~~~~~~~~~~~~  Create a bunch of individual debris footprints ~~~~~~~~~~~~~~~~~~~~~~~
//        // Generate a bunch of footprints at different explosion times
//        string CapeGoogleEarthFileITER("GeneratedFiles/CapeFootPrint");
//        string datExtension(".dat");
//        string kmlExtension(".kml");
//        
//        
//        for (double timeExplodeUntil = 70; timeExplodeUntil < 71.; timeExplodeUntil += 10.){
//            char datfname[CapeGoogleEarthFileITER.size() + 2 + 4];
//            sprintf(datfname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, datExtension.c_str() );
//            
//            char kmlfname[CapeGoogleEarthFileITER.size() + 2 + 4];
//            sprintf(kmlfname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, kmlExtension.c_str() );
//            
//            //        sprintf(fname, "%s%d%s",  altExplodeUntil, extension.c_str() );
//            cout << "datfname = " << datfname << endl;
//            cout << "kmlfname = " << kmlfname << endl;
//            
//            {
//                char StateOption[] = "full rocket";
//                Trajectory MyCape(StateOption, WindOption, DensityOption);
//                MyCape.Initialize_Stages(numPerBatch, delta_t, cape_nominal_traj_filename);  //seed is still fixed!!!
//                MyCape.PropagateWholeRocket();
//                int timesteupsUntilCutoff = MyCape.getTimestepsUntilCutoff();
//                MyCape.write_to_google_earth_native(CapeGoogleEarthTrajFile, numToWrite);   // Makes a google earth movie of TRAJECTORIES, not compact envelope
//                MyCape.write_LRHC_points_file(CapeLrhcFile, tstepMinutes);              // Exports rocket launch points as lat lon for making footprint
//                MyCape.write_stage_down_points_file(CapeStageDownFile, tstepMinutes);
//                
//                Footprint3D MyCapeFootprint(CapeLrhcFile, binSizeKm);
//                MyCapeFootprint.RemoveFootprintKeepTiming();
//                MyCape.MonteCarloDebris(MyCapeFootprint, timeExplodeUntil);
//                cout << "storing as vector~~~~~~~~~~~~~~" << endl;
//                MyCapeFootprint.store_footprint_as_vector(datfname);
//                //                cout << "exporting~~~~~~~~~~~~~~~~~~~~" << endl;
//                //                MyCapeFootprint.exportGoogleEarth(kmlfname);
//                
//            } }
//        
//        // ~~~~~~~~~~~~~~~~ End individual debris clouds genreation ~~~~~~~~~~~~~~~~~~~~~

        
        
        
        
        
        // Do a launch from the Cape
//        char StateOption[] = "full rocket";
//        Trajectory MyCape(StateOption, WindOption, DensityOption);
//        MyCape.Initialize_Stages(numPerBatch, delta_t, cape_nominal_traj_filename);  //seed is still fixed!!!
//        MyCape.PropagateWholeRocket();
//        int timesteupsUntilCutoff = MyCape.getTimestepsUntilCutoff();
//        MyCape.write_to_google_earth_native(CapeGoogleEarthTrajFile, numToWrite);   // Makes a google earth movie of TRAJECTORIES, not compact envelope
//        MyCape.write_LRHC_points_file(CapeLrhcFile, tstepMinutes);              // Exports rocket launch points as lat lon for making footprint
//        MyCape.write_stage_down_points_file(CapeStageDownFile, tstepMinutes);   // Exports falling first stage as points for making footprint

        // Basic footprint stuff
//        Footprint3D MyCapeFootprint(CapeLrhcFile, binSizeKm);
//        MyCape.MonteCarloDebris(MyCapeFootprint, 10.);
//        MyCapeFootprint.exportGoogleEarth(CapeGoogleEarthFile);
//        //        MyCape.HowLongUntilEnterNAS();
        
        // Basic skygrid stuff
//        SkyGrid CapeSkyGrid(CapeLrhcFile, xBinLength, yBinLength, zBinHeight);
//        MyCape.MonteCarloDebris(CapeSkyGrid);
//        MyCape.HowLongUntilEnterNAS();
//        CapeSkyGrid.getKDE(0.5,3);
        
        
        
//        // ~~~~~~~~~~~~~~~~ Basics of running the SkyGrid KDE ~~~~~~~~~~~~~~~~~~~~~~~~~~
//        
//        char SkyGridVectorFile[] = "GeneratedFiles/SkyGridVector.dat";
//        char TestGEFootprint[] = "GeneratedFiles/TestGoogEarthFP.kml";
//
////        // Do a launch from the Cape
////        char StateOption[] = "full rocket";
////        Trajectory MyCape(StateOption, WindOption, DensityOption);
////        MyCape.Initialize_Stages(numPerBatch, delta_t, cape_nominal_traj_filename);  //seed is still fixed!!!
////        MyCape.PropagateWholeRocket();
////        int timesteupsUntilCutoff = MyCape.getTimestepsUntilCutoff();
////        MyCape.write_to_google_earth_native(CapeGoogleEarthTrajFile, numToWrite);   // Makes a google earth movie of TRAJECTORIES, not compact envelope
////        MyCape.write_LRHC_points_file(CapeLrhcFile, tstepMinutes);              // Exports rocket launch points as lat lon for making footprint
////        MyCape.write_stage_down_points_file(CapeStageDownFile, tstepMinutes);   // Exports falling first stage as points for making footprint
////        
////        // Basic skygrid stuff
////        double timeExplodUntil = 60.;
////        SkyGrid CapeSkyGrid(CapeLrhcFile, xBinLength, yBinLength, zBinHeight);
////        MyCape.MonteCarloDebris(CapeSkyGrid, timeExplodUntil);
////        MyCape.HowLongUntilEnterNAS();
////        
////        CapeSkyGrid.StoreGrid(SkyGridVectorFile);
//        
//        double deltaXY = 0.5;
//        SkyGrid NewSkyGrid(SkyGridVectorFile);
//        NewSkyGrid.generateAllPointsFromKDE(deltaXY);
//        Footprint3D FpFromSkyGrid(NewSkyGrid);
//        FpFromSkyGrid.exportGoogleEarth(TestGEFootprint);
//        
//        // ~~~~~~~~~~~~~~~ End of SkyGrid KDE Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
//        // ~~~~~~~~~~~~~~~  Create a bunch of individual debris footprints WITH SKYGRID ~~~~~~~~~~~~~~~~~~~~~~~
//        // Generate a bunch of footprints at different explosion times
//        char StateOption[] = "full rocket";
//        string CapeGoogleEarthFileITER("GeneratedFiles/CapeFootPrint");
////        string CapeGoogleEarthFileITER("GeneratedFiles/CapeFootPrint");
////        string extension(".dat");
//        string extension(".kml");
//
//
//        double deltaXY = 0.5;
//        
//        for (double timeExplodeUntil = 150; timeExplodeUntil < 151.; timeExplodeUntil += 10.){
//            char fname[CapeGoogleEarthFileITER.size() + 2 + 4];
//            sprintf(fname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, extension.c_str() );
//            //        sprintf(fname, "%s%d%s",  altExplodeUntil, extension.c_str() );
//            cout << "fname = " << fname << endl;
//
//            {
//                Trajectory MyCape(StateOption, WindOption, DensityOption);
//                MyCape.Initialize_Stages(numPerBatch, delta_t, cape_nominal_traj_filename);  //seed is still fixed!!!
//                MyCape.PropagateWholeRocket();
//                int timesteupsUntilCutoff = MyCape.getTimestepsUntilCutoff();
////                MyCape.write_to_google_earth_native(CapeGoogleEarthTrajFile, numToWrite);   // Makes a google earth movie of TRAJECTORIES, not compact envelope
//                MyCape.write_LRHC_points_file(CapeLrhcFile, tstepMinutes);              // Exports rocket launch points as lat lon for making footprint
//                MyCape.write_stage_down_points_file(CapeStageDownFile, tstepMinutes);
//
////                Footprint3D MyCapeFootprint(CapeLrhcFile, binSizeKm);
////                MyCapeFootprint.RemoveFootprintKeepTiming();
////                MyCape.MonteCarloDebris(MyCapeFootprint, timeExplodeUntil);
////                MyCapeFootprint.store_footprint_as_vector(fname);
//                
////                SkyGrid NewSkyGrid(SkyGridVectorFile);
//                SkyGrid CapeSkyGrid(CapeLrhcFile, xBinLength, yBinLength, zBinHeight);
//                CapeSkyGrid.RemoveHistogramKeepTiming();
//                MyCape.MonteCarloDebris(CapeSkyGrid, timeExplodeUntil);
//                CapeSkyGrid.generateAllPointsFromKDE(deltaXY);
//                Footprint3D FpFromSkyGrid(CapeSkyGrid);
//                FpFromSkyGrid.exportGoogleEarth(fname);
//
//
//            }
//        }
//        
//        // ~~~~~~~~~~~~~~~~ End individual debris clouds genreation WITH SKYGRID ~~~~~~~~~~~~~~~~~~~~~
        
        
        
        
        
        

        
        
        
        
        
        
        
        
//        // ~~~~~~~~~~~~~  Generate Suborbital Footprint ~~~~~~~~~~~~~~~~~~~~~
//        // Set Atmosphere Options
//        char StateOption[] = "suborbital";
//        string nominal_traj_filename = "Files/SS2/SS2Acceleration.txt";
//        double delta_t = 2.0;
//
//        char Ss2LrhcFile[] = "GeneratedFiles/Ss2LRHC.dat";
//        char Ss2GoogleEarth[] = "GeneratedFiles/Ss2Footprint.kml";
//        
//        Trajectory MyTraj(StateOption, WindOption, DensityOption);
//        MyTraj.Initialize_Suborbital(delta_t, nominal_traj_filename);
//        MyTraj.Propagate_Suborbital();
////        MyTraj.write_to_google_earth_native(Ss2GoogleEarthTrajFile, numToWrite);
//        MyTraj.write_LRHC_points_file(Ss2LrhcFile, tstepMinutes);              // Exports rocket launch points as lat lon for making footprint
//        Footprint3D Ss2Footprint(Ss2LrhcFile, binSizeKm);
//        MyTraj.MonteCarloDebris(Ss2Footprint, 60.);
//
//        Ss2Footprint.exportGoogleEarth(Ss2GoogleEarth);
//        // ~~~~~~~~~~~~~~ Done With Suborbital ~~~~~~~~~~~~~~~~~~~~~~~~~~~


        
    

        
//    } else {
    

//        // ~~~~~~~~~~~ Create an updating footprint ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//        Footprint3D updatingFootprintZeroToOne("GeneratedFiles/CapeFootPrint 10.dat");
//
//        // Generate a bunch of footprints at different explosion times
//        string CapeGoogleEarthFileITER("GeneratedFiles/CapeFootPrint");
//        string extension(".dat");
//
//        for (double timeExplodeUntil = 20; timeExplodeUntil < 131.; timeExplodeUntil += 10.){
//            char fname[CapeGoogleEarthFileITER.size() + 3 + 4];
//            sprintf(fname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, extension.c_str() );
//            cout << "fname = " << fname << endl;
//
//            Footprint3D tempFootprint(fname);
//            updatingFootprintZeroToOne.MergeFootprintVectors(tempFootprint);
//
//        }
//        cout << "launcLatTEMP = " << updatingFootprintZeroToOne.getLaunchLat() << endl;
//        
//        updatingFootprintZeroToOne.ProjectAllPointsDown();
//        
//        cout << "launcLatTEMP = " << updatingFootprintZeroToOne.getLaunchLat() << endl;
//        
//        updatingFootprintZeroToOne.ChopTimeAt(1);
//
//
//        Footprint3D updatingFootprintOneToTwo("GeneratedFiles/CapeFootPrint 60.dat");
//        for (double timeExplodeUntil = 70; timeExplodeUntil < 131.; timeExplodeUntil += 10.){
//            char fname[CapeGoogleEarthFileITER.size() + 3 + 4];
//            sprintf(fname, "%s%3.0f%s", CapeGoogleEarthFileITER.c_str(), timeExplodeUntil, extension.c_str() );
//            cout << "fname = " << fname << endl;
//
//            Footprint3D tempFootprint(fname);
//            updatingFootprintOneToTwo.MergeFootprintVectors(tempFootprint);
//
//        }
//        updatingFootprintOneToTwo.ProjectAllPointsDown();
//        updatingFootprintOneToTwo.ChopTimeAt(1);    //makes this fp last for one minute
//        updatingFootprintOneToTwo.AddToFootprintUTC(0, 1);  //plus hours,minutes
//
//        Footprint3D updatingFootprint;
//        updatingFootprint = updatingFootprintZeroToOne;
//        updatingFootprint.MergeFootprintVectors(updatingFootprintOneToTwo);
//        
//        // Translate footprint to Georgia  30.921463,-81.513451
//        updatingFootprint.ChangeLaunchSiteToDeg(30.921463, -81.513451);
//        updatingFootprint.exportGoogleEarth(CapeGoogleEarthFile);
//
//        char CapeGoogleEarthFile2[] = "GeneratedFiles/CapeFootPrintBlowStaging.kml";
//        Footprint3D FootprintStaging("GeneratedFiles/CapeFootPrint170.dat");
//        FootprintStaging.ChangeLaunchSiteToDeg(30.921463, -81.513451);
//        FootprintStaging.exportGoogleEarth(CapeGoogleEarthFile2);
//
//
//        char CapeGoogleEarthFile3[] = "GeneratedFiles/CapeFootPrintBlowMaxQ.kml";
//        Footprint3D FootprintMaxQ("GeneratedFiles/CapeFootPrint 90.dat");
//        FootprintMaxQ.exportGoogleEarth(CapeGoogleEarthFile3);
//
//        char CapeGoogleEarthFile4[] = "GeneratedFiles/CapeFootPrintLRHC.kml";
//        Footprint3D FootprintLHRC(CapeLrhcFile, binSizeKm);
//        Footprint3D FootprintEarlyDebris("GeneratedFiles/CapeFootPrint 70.dat");
//        FootprintLHRC.MergeFootprintVectors(FootprintEarlyDebris);
//        FootprintLHRC.exportGoogleEarth(CapeGoogleEarthFile4);
//        
//        // ~~~~~~~~~~~~~~ End create updating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
        
        
        
        
        
        // The trajectory footprint
//        Footprint3D CapeFootprint(CapeFootprintVectorFile);
        
//        CapeFootprint.ChangeLaunchSiteTo(37.8338*(PI/180), -75.4882*(PI/180));
//        CapeFootprint.ChangeLaunchSiteToDeg(32.8894, -106.9995);
//        CapeFootprint.ChangeAzimuthToDeg(90);
        
//        cout << "current footprint launchLat = " << CapeFootprint.getLaunchLat() << endl;

//        CapeFootprint.ShiftFootprintByMinutes(40);
        
        
        
        
        
//        // First stretch out Cape
//        Footprint3D shiftedFootprint;
//        shiftedFootprint = CapeFootprint;
//        
//        int shiftHowManyMinutes = 40;
//        for (int minCount = 0; minCount < shiftHowManyMinutes; minCount++){
//            shiftedFootprint.SetFootprintUTC(0, 1);
//            CapeFootprint.MergeFootprintVectors(shiftedFootprint); }
//
//        // Combine them all into something smooth
//        CapeFootprint.SmoothedOut();
        
        
        
//        // Put in the fallen stage
//        Footprint3D MyStageDown(CapeStageDownFile, binSizeKm);
//        Footprint3D shiftedStageDown;
//        shiftedStageDown = MyStageDown;
//        
//        for (int minCount = 0; minCount < shiftHowManyMinutes; minCount++){
//            shiftedStageDown.SetFootprintUTC(0, 1);
//            MyStageDown.MergeFootprintVectors(shiftedStageDown); }
//        
//        MyStageDown.SmoothedOut();
//        
//        // Put em all together
//        CapeFootprint.MergeFootprintVectors(MyStageDown);
        
        
//        // Now because Google Earth sucks a little at movies, simulate a loop myself
//        shiftedFootprint = CapeFootprint;
//        shiftHowManyMinutes = 6;
//        for (int minCount = 0; minCount < shiftHowManyMinutes; minCount++){
//            shiftedFootprint.SetFootprintUTC(0, 90);
//            CapeFootprint.MergeFootprintVectors(shiftedFootprint); }
        
        
        
//        Footprint3D MarsFootprint(MarsFootprintVectorFile);
//        double curUTC = CapeFootprint.GetFootprintUTC();
//        MarsFootprint.SetFootprintUTC(curUTC);
        
//        MarsFootprint.ProjectAllPointsDown();
//
        // Make it launch 45 minutes later
        
        //105 timesteps
//        MarsFootprint.ChopTimeAt(3.5);  //minutes
        
        
        
        
        
        
        
//        MarsFootprint.AddToFootprintUTC(0, 45);
//        MarsFootprint.ShiftFootprintByMinutes(45);

//
//        // Now stretch out Mars
//        shiftedFootprint = MarsFootprint;
//        
//        shiftHowManyMinutes = 20;
//        for (int minCount = 0; minCount < shiftHowManyMinutes; minCount++){
//            shiftedFootprint.SetFootprintUTC(0, 1);
//            MarsFootprint.MergeFootprintVectors(shiftedFootprint); }
//        
//        // Combine them all into something smooth
//        MarsFootprint.SmoothedOut();
        
        
        // Now put them together
//        Footprint3D TotalMission;
//        TotalMission = CapeFootprint;
//        TotalMission.MergeFootprintVectors(MarsFootprint);

//        Footprint3D TotalMission;
////        TotalMission = updatingFootprint;
//        TotalMission = FootprintStaging;
//        
//        
////        TotalMission = MarsFootprint;
//        
//        
////        MyFootprint.store_footprint_as_points();
////        TotalMission.exportGoogleEarth(FullMissionGoogleEarthFile);
//
//        
//        // Make FACET files
//        bool makeFacet = true;
//        if (makeFacet) {
//            int launchTimeHours = 9;
//            int launchTimeMinutes = 45;
//            int offsetTimeMinutes = 0;	//This turns on the first SUA $offset minutes earlier than the launch time
//            
//            
//            string folderPath("GeneratedFacet/");
//            string folderBaseName("Sandbox");
//            char buffer[40];
//            sprintf(buffer,"__Time_%02d%02d", (int) floor(launchTimeHours),(int) floor(launchTimeMinutes));
//            string timeStartInfo(buffer);
//
//            sprintf(buffer, "__dTime_%02d", (int) floor(tstepMinutes));
//            string timeDeltaInfo(buffer);
//
//            sprintf(buffer,"__BinSizeKm_%.2f", binSizeKm);
//            string binInfo(buffer);
//
//            string folderName = folderPath + folderBaseName + timeStartInfo + timeDeltaInfo + binInfo + "/";
//
//            string mkdir("mkdir ");
//            system((mkdir + folderName).c_str());
//            
//            cout << (mkdir + folderName).c_str() << endl;
//
//            TotalMission.make_facet_files(folderName, launchTimeHours*60 + launchTimeMinutes, offsetTimeMinutes, tstepMinutes);
//        }
    
        
//    }
    

//	string endMessage("say \"finished\"");
//	system(endMessage.c_str());
//
//    return 0;
//}


//	cout << "NOTE: Staging is not really handled properly at the moment, so don't trust anything beyond the first stage" << endl;
//
//	// Actually, it probably doesn't matter that much because by the time we're making big latitude changes, we'll be beyond the regime i care about.
//	cout << "NOTE: Calculations of altitude-dependent quantities is not quite right.  Need to properly implement oblate earth radius" << endl;
//	cout << "NOTE: Passing a 'Debris' object but didn't write a copy constructor: BE WARY OF ERRORS HERE" << endl;
//
//	cout << "NOTE: Need to do something about availablepointstorage and whatnot" << endl;













//int main (int argc, char * const argv[]) {
//	
//	int WhatToDo = Do_Propagate;
//    //	int WhatToDo = Do_SimLeftRightHotCold;
//    //	int WhatToDo = Do_MakeFacetFiles;
//    //	int WhatToDo = Do_MC;
//    //	int WhatToDo = Merge_Footprint_Vectors;
//	
//	double tstepMinutes = 1.;
//	double binSizeKm = 1.;
//	int launchTimeHours = 15.;
//	int launchTimeMinutes = 45.;
//	int offsetTimeMinutes = 5;	//This turns on the first SUA $offset minutes earlier than the launch time
//	
//    //	double debrisCutoff = 20.;	//[km]
//	
//	string nominalFileName("GeneratedFiles/trajectory_points.dat");
//	string rawPointsName("GeneratedFiles/points_to_wrap_up_abridged.dat");
//	string lrhcFileName("GeneratedFiles/points_to_wrap_up_LRHC.dat");
//	
//	// Set Atmosphere Options
//	char WindOption[] = "simple atmosphere";
//	char DensityOption[] = "cantwell density";
//    
//	char StateOption[] = "first stage";
//	string nominal_traj_filename = "Files/zzz_nominal_cape_finished.txt";
//    
//    //    char StateOption[] = "suborbital";
//    //	string nominal_traj_filename = "Files/SS2/SS2Acceleration.txt";
//    
//	string footprintVectorOutputFile = "GeneratedFiles/CAPE_Footprint_lowBallistic.dat";
//    
//    
//    //	string nominal_traj_filename = "Files/zzz_nominal_MARS.txt";
//    //	string footprintVectorOutputFile = "GeneratedFiles/MARS_Footprint_lowBallistic.dat";
//    exit(90);
//    
//    
//	switch (WhatToDo) {
//            
//            //		case Do_SimLeftRightHotCold:{
//            //			Architecture LeftRightHotCold;
//            //			int fidelity = 30;
//            //			double semiMajorPercent = 2.;
//            //			double eccentricity = 0.95;
//            //			LeftRightHotCold.SimLeftRightHotCold(nominalFileName.c_str(), fidelity, semiMajorPercent, eccentricity);
//            //			// USE write_LRHC_points
//            ////			LeftRightHotCold.write_all_points(lrhcFileName);
//            //			}break;
//            
//            
//		case Do_MC:{
//			// Propagate Debris Assuming Explosion At Every Point
//			
//			timer clock2;
//			clock2.start();
//			Architecture SimpleArch;
//			SimpleArch.MonteCarlo(WindOption, DensityOption, nominal_traj_filename, footprintVectorOutputFile);
//			clock2.stop();
//            
//			cout << "That took " << clock2.how_long() << " seconds" << endl;
//			
//        }break;
//			
//		case Merge_Footprint_Vectors:{
//			string FootprintVector2 = "GeneratedFiles/MARS_Footprint_lowBallistic.dat";
//			string FootprintVector1 = "GeneratedFiles/CAPE_Footprint_lowBallistic.dat";
//			string FootprintVectorOutput = "GeneratedFiles/Cape_Mars_Footprint.dat";
//			
//			cout << "Goign to merge that footprint~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
//			Footprint3D MyFootprint(FootprintVector1);
//			MyFootprint.SetFootprintUTC(0.,30.);	//add some hours and minutes
//			MyFootprint.MergeFootprintVectors(FootprintVector2, FootprintVectorOutput);
//			// Now take that new footprint, translate it forward, and merge it with its old self
//            //			MyFootprint.SetFootprintUTC(1.,0);	//add some hours and minutes
//            //			MyFootprint.MergeFootprintVectors(FootprintVectorOutput, FootprintVectorOutput);
//            
//            
//			
//			cout << "done merging" << endl;
//            
//            //			string folderPath("/Users/marian/Documents/Research/GeneratedFiles/");
//            //			MyFootprint.make_facet_files(folderName, launchTimeHours*60 + launchTimeMinutes, offsetTimeMinutes, tstepMinutes);
//			// DEBUGGING
//			double launchTimeHours = 15.;
//			double launchTimeMinutes = 45.;
//			double offsetTimeMinutes = 0.;
//            //			string folderName("/Users/marian/Documents/Research/GeneratedFiles/DEBUG/");
//			string folderName("/Users/marian/Dropbox/To_Facet/DEBUG/");
//            
//			MyFootprint.make_facet_files(folderName, launchTimeHours*60 + launchTimeMinutes, offsetTimeMinutes, tstepMinutes);
//			// END DEBUGGING
//			
//        }break;
//			
//        case Do_MakeFacetFiles:{
//            string pointsFile("GeneratedFiles/points_to_wrap_up_LRHC.dat");
//            //				string pointsFile("GeneratedFiles/points_to_wrap_up_abridged.dat");
//            string folderPath("/Users/marian/Documents/Research/GeneratedFiles/");
//            //			string folderBaseName("All\\ Points\\ Explode\\ cutoff20");
//            string folderBaseName("Sandbox");
//            char buffer[40];
//            sprintf(buffer,",\\ Time\\ %02d%02d", (int) floor(launchTimeHours),(int) floor(launchTimeMinutes));
//            string timeStartInfo(buffer);
//            
//            sprintf(buffer, ",\\ dTime\\ %2d", (int) floor(tstepMinutes));
//            string timeDeltaInfo(buffer);
//            
//            sprintf(buffer,",\\ BinSizeKm\\ %.2f", binSizeKm);
//            string binInfo(buffer);
//            
//            string folderName = folderPath + folderBaseName + timeStartInfo + timeDeltaInfo + binInfo + "/";
//            
//            string mkdir("mkdir ");
//            system((mkdir + folderName).c_str());
//            Footprint3D my_footprint(pointsFile,binSizeKm);
//            my_footprint.generate_footprint_at_timesteps();
//            my_footprint.make_facet_files(folderName, launchTimeHours*60 + launchTimeMinutes, offsetTimeMinutes, tstepMinutes);
//            
//            string googleEarthFile("Envelope.kml");
//            my_footprint.exportGoogleEarth(googleEarthFile);
//        }break;
//            
//        case Do_Propagate:{
//            // Set Misc Options
//            int num_launches_per_batch = 1;
//            double delta_t = 2.;
//            
//            Trajectory MyTraj(StateOption, WindOption, DensityOption);
//            MyTraj.Initialize_First_Stage(num_launches_per_batch, delta_t, nominal_traj_filename);
//            MyTraj.Propagate_First_Stage();
//            
//            //				MyTraj.Initialize_Suborbital(delta_t, nominal_traj_filename);
//            //                MyTraj.Propagate_Suborbital();
//            
//            char outFileName[] = "Files/cppOutFile.dat";
//            //MyTraj.write_batch_to_google_earth("useless", 0);     //CANDIDATE FOR DELETION
//            MyTraj.write_single_trajectory_to_file(outFileName);	//Should specify the file name!!!
//            MyTraj.write_to_google_earth_native("useless",666);
//            // Propagate Debris Assuming Explosion At Every Point
//            cout << "Exiting early" << endl;
//            exit(8);
//            //				Architecture SimpleArch;
//            //				SimpleArch.read_single_trajectory_from_file(nominalFileName);	//Should specify the file name!!!
//            //				SimpleArch.explode_at_all_points();			//Should specify the atmosphere options!!!
//            
//            
//            
//            
//            
//            
//        }break;
//            
//            
//		default:{
//			cout << "ERROR!  You chose some WhatToDo option that doesn't exist ~~~~~~~~~~~~~~~~~~" << endl;
//        }break;
//	}
//	
//    
//	
//	cout << "NOTE: Staging is not really handled properly at the moment, so don't trust anything beyond the first stage" << endl;
//	
//	// Actually, it probably doesn't matter that much because by the time we're making big latitude changes, we'll be beyond the regime i care about.
//	cout << "NOTE: Calculations of altitude-dependent quantities is not quite right.  Need to properly implement oblate earth radius" << endl;
//	cout << "NOTE: Passing a 'Debris' object but didn't write a copy constructor: BE WARY OF ERRORS HERE" << endl;
//    
//	cout << "NOTE: Need to do something about availablepointstorage and whatnot" << endl;
//	
//    
//	string endMessage("say \"finished\"");
//	system(endMessage.c_str());
//    
//	
//	
//	
//	
//    return 0;
//}











































