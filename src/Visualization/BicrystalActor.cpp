/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BicrystalActor_cpp_
#define model_BicrystalActor_cpp_

#include <iostream>
#include <deque>
#include <string>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>

#include <TerminalColors.h>
#include <BicrystalActor.h>

namespace gbLAB
{
    
        
        /**********************************************************************/
        BicrystalActor::BicrystalActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const rndr) :
        /* init */ renderWindow(renWin)
        /* init */,renderer(rndr)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showA(new QCheckBox(this))
        /* init */,showB(new QCheckBox(this))

//        /* init */,showNodeLabels(new QCheckBox(this))
//        /* init */,showVelocities(new QCheckBox(this))
//        /* init */,velocityScaleEdit(new QLineEdit("1"))
        /* init */,aPolyData(vtkSmartPointer<vtkPolyData>::New())
        /* init */,aGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
        /* init */,aMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
        /* init */,aActor(vtkSmartPointer<vtkActor>::New())
/* init */,bPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,bGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
/* init */,bMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,bActor(vtkSmartPointer<vtkActor>::New())

//        /* init */,labelPolyData(vtkSmartPointer<vtkPolyData>::New())
//        /* init */,labelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New())
//        /* init */,labelActor(vtkSmartPointer<vtkActor2D>::New())
//        /* init */,velocityPolyData(vtkSmartPointer<vtkPolyData>::New())
//        /* init */,velocityGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
//        /* init */,velocityMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
//        /* init */,velocityActor(vtkSmartPointer<vtkActor>::New())
//        /* init */,singleNodeLabelPolyData(vtkSmartPointer<vtkPolyData>::New())
//        /* init */,singleNodeLabelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New())
//        /* init */,singleNodeLabelActor(vtkSmartPointer<vtkActor2D>::New())
//        /* init */,singleNodeID(0)
//        /* init */,nodeClr{{100,100,100},{0,255,255},{255,0,255},{1,1,1}}
        {
            
            showA->setText("lattice A");
            showA->setChecked(true);
            aActor->SetVisibility(true);
            showB->setText("lattice B");
            showB->setChecked(true);
            bActor->SetVisibility(true);

//            showNodeLabels->setText("node labels");
//            showNodeLabels->setChecked(false);
//            labelActor->SetVisibility(false);
//
//            showVelocities->setText("velocities");
//            showVelocities->setChecked(false);
//            velocityActor->SetVisibility(false);
//            velocityScaleEdit->setEnabled(false);


            mainLayout->addWidget(showA,0,0,1,1);
            mainLayout->addWidget(showB,1,0,1,1);

//            mainLayout->addWidget(showNodeLabels,1,0,1,1);
//            mainLayout->addWidget(showVelocities,2,0,1,1);
//            mainLayout->addWidget(velocityScaleEdit,2,1,1,1);
            this->setLayout(mainLayout);

            connect(showA,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(showB,SIGNAL(stateChanged(int)), this, SLOT(modify()));

//            connect(showNodeLabels,SIGNAL(stateChanged(int)), this, SLOT(modify()));
//            connect(showVelocities,SIGNAL(stateChanged(int)), this, SLOT(modify()));
//            connect(velocityScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));

            
            aGlyphs->SetInputData(aPolyData);
            aGlyphs->SetSourceConnection(vtkSmartPointer<vtkSphereSource>::New()->GetOutputPort());
            aGlyphs->ScalingOn();
            aGlyphs->SetScaleModeToScaleByVector();
            aGlyphs->SetScaleFactor(0.25);
            aGlyphs->SetColorModeToColorByScalar();
            aGlyphs->Update();
            aMapper->SetInputConnection(aGlyphs->GetOutputPort());
            aActor->SetMapper(aMapper);
            aActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
            renderer->AddActor(aActor);
            
            bGlyphs->SetInputData(bPolyData);
            bGlyphs->SetSourceConnection(vtkSmartPointer<vtkSphereSource>::New()->GetOutputPort());
            bGlyphs->ScalingOn();
            bGlyphs->SetScaleModeToScaleByVector();
            bGlyphs->SetScaleFactor(0.25);
            bGlyphs->SetColorModeToColorByScalar();
            bGlyphs->Update();
            bMapper->SetInputConnection(bGlyphs->GetOutputPort());
            bActor->SetMapper(bMapper);
            bActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
            renderer->AddActor(bActor);

            
            // Labels
//            labelMapper->SetInputData(labelPolyData);
//            labelMapper->SetLabelModeToLabelScalars();
//            labelMapper->SetLabelFormat("%1.0f");
//            labelMapper->GetLabelTextProperty()->SetFontSize(20);
//            labelActor->SetMapper(labelMapper);
//            labelActor->GetProperty()->SetColor(0.0, 0.0, 0.0); //(R,G,B)
//
//            // Velocities
//            velocityGlyphs->SetInputData(velocityPolyData);
//            velocityGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
//            velocityGlyphs->ScalingOn();
//            velocityGlyphs->SetScaleModeToScaleByVector();
//            velocityGlyphs->OrientOn();
//            velocityGlyphs->ClampingOff();
//            velocityGlyphs->SetVectorModeToUseVector();
//            velocityGlyphs->SetIndexModeToOff();
//            velocityMapper->SetInputConnection(velocityGlyphs->GetOutputPort());
//            velocityMapper->ScalarVisibilityOff();
//            velocityActor->SetMapper(velocityMapper);
//            velocityActor->GetProperty()->SetColor(1.0, 0.0, 1.0); //(R,G,B)
//
//            // Single node Label
//            singleNodeLabelMapper->SetInputData(singleNodeLabelPolyData);
//            singleNodeLabelMapper->SetLabelModeToLabelScalars();
//            singleNodeLabelMapper->SetLabelFormat("%1.0f");
//            singleNodeLabelMapper->GetLabelTextProperty()->SetFontSize(20);
//            singleNodeLabelActor->SetMapper(singleNodeLabelMapper);
//            singleNodeLabelActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
//            singleNodeLabelActor->VisibilityOff();
            
//            renderer->AddActor(velocityActor);
//            renderer->AddActor(labelActor);
//            renderer->AddActor(singleNodeLabelActor);

        }
        

        
        /**********************************************************************/
        void BicrystalActor::updateConfiguration(const std::shared_ptr<BiCrystal<3>>& bc)
        {// https://stackoverflow.com/questions/6878263/remove-individual-points-from-vtkpoints
            std::cout<<"Updating BiCrystal..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
//
            vtkSmartPointer<vtkPoints> aPoints(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkPoints> bPoints(vtkSmartPointer<vtkPoints>::New());
            
            const int N=5;
            for(int i=-N;i<N+1;++i)
            {
                for(int j=-N;j<N+1;++j)
                {
                    for(int k=-N;k<N+1;++k)
                    {
                        const Eigen::Matrix<double,3,1> Pa(bc->A.latticeBasis*(Eigen::Matrix<double,3,1>()<<i,j,k).finished());
                        aPoints->InsertNextPoint(Pa.data());
                        const Eigen::Matrix<double,3,1> Pb(bc->B.latticeBasis*(Eigen::Matrix<double,3,1>()<<i,j,k).finished());
                        bPoints->InsertNextPoint(Pb.data());

                    }
                }
            }
            

//            vtkSmartPointer<vtkUnsignedCharArray> nodeColors(vtkSmartPointer<vtkUnsignedCharArray>::New());
//            nodeColors->SetNumberOfComponents(3);
//
//            vtkSmartPointer<vtkDoubleArray> nodeLabels(vtkSmartPointer<vtkDoubleArray>::New());
//            nodeLabels->SetNumberOfComponents(1);
//            nodeLabels->SetName("node IDs");
//
//            vtkSmartPointer<vtkPoints> singleNodePoint(vtkSmartPointer<vtkPoints>::New());
//            vtkSmartPointer<vtkDoubleArray> singleNodenodeLabels(vtkSmartPointer<vtkDoubleArray>::New());
//
//            vtkSmartPointer<vtkDoubleArray> velocityVectors(vtkSmartPointer<vtkDoubleArray>::New());
//            velocityVectors->SetNumberOfComponents(3);
//            velocityVectors->SetName("nodeVelocity");
//
//            for(const auto& node : configIO.nodes())
//            {
//                nodePoints->InsertNextPoint(node.P.data());
//                nodeLabels->InsertNextTuple1(node.sID);
//                velocityVectors->InsertNextTuple(node.V.data()); // arrow vector
//                nodeColors->InsertNextTypedTuple(node.meshLocation>2? this->nodeClr[3] : this->nodeClr[node.meshLocation]);
//
//                // Single node
//                if(node.sID==singleNodeID)
//                {
//                    singleNodePoint->InsertNextPoint(node.P.data());
//                    singleNodenodeLabels->InsertNextTuple1(node.sID);
//                }
//            }
//
            aPolyData->SetPoints(aPoints);
            aPolyData->Modified();
            bPolyData->SetPoints(bPoints);
            bPolyData->Modified();

//
//            labelPolyData->SetPoints(nodePoints);
//            labelPolyData->GetPointData()->SetScalars(nodeLabels);
//            labelPolyData->Modified();
//
//            singleNodeLabelPolyData->SetPoints(singleNodePoint);
//            singleNodeLabelPolyData->GetPointData()->SetScalars(singleNodenodeLabels);
//            singleNodeLabelPolyData->Modified();
//
//            velocityPolyData->SetPoints(nodePoints);
//            velocityPolyData->GetPointData()->SetVectors(velocityVectors);
//            velocityPolyData->Modified();
            renderer->ResetCamera();
            renderWindow->Render();
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void BicrystalActor::modify()
        {
            
            aActor->SetVisibility(showA->isChecked());
            bActor->SetVisibility(showB->isChecked());

//            labelActor->SetVisibility(showNodeLabels->isChecked());
//            velocityActor->SetVisibility(showVelocities->isChecked());
//            velocityScaleEdit->setEnabled(showVelocities->isChecked());
//            const double vScaling(std::atof(velocityScaleEdit->text() .toStdString().c_str()));
//            velocityGlyphs->SetScaleFactor(vScaling);

//                velocityActor->SetScale(vScaling);
//            aGlyphs->SetScaleFactor(2.0*this->tubeRadius*1.2);
//            
//            if(this->showVelocities)
//            {
//                velocityActor->VisibilityOn();
//            }
//            else
//            {
//                velocityActor->VisibilityOff();
//            }
//            
//            if(this->showNodeIDs)
//            {
//                labelActor->VisibilityOn();
//                
//            }
//            else
//            {
//                labelActor->VisibilityOff();
//            }
//            
//            velocityGlyphs->SetScaleFactor(this->velocityFactor);
//            
//            
//            if(this->showSingleNode)
//            {
//                // HERE WE SHOULD CHANGE THE NODE POSITION BASED ON NODE ID
//                // OTHERWISE THE SELECTED NODE WILL BE VISIBLE ONLY UPON LOADING A NEW FRAME
//                std::cout<<"RELOAD FRAME TO SHOW SELECTED NODE"<<std::endl;
//                singleNodeLabelActor->VisibilityOn();
//            }
//            else
//            {
//                singleNodeLabelActor->VisibilityOff();
//            }
//            
//            if(this->showA)
//            {
//                aActor->VisibilityOn();
//            }
//            else
//            {
//                aActor->VisibilityOff();
//            }

            renderWindow->Render();
        }
        
    
} // namespace model
#endif
