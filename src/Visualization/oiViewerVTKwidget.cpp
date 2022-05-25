/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_oiViewerVTKwidget_cpp_
#define model_oiViewerVTKwidget_cpp_

#include <QFileDialog>

#include <oiViewerVTKwidget.h>


namespace gbLAB
{
    


    
    oiViewerVTKwidget::oiViewerVTKwidget(QWidget*) :
    /* init */ mainLayout(new QGridLayout(this))
    /* init */,loadBicrystalButton(new QPushButton("Load",this))
    /* init */,tabWidget(new QTabWidget(this))
    /* init */,openglWidget(new QVTKOpenGLStereoWidget(this))
    /* init */,renderWindow(vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New())
    /* init */,renderer(vtkSmartPointer<vtkRenderer>::New())
//    /* init */,traitsIO(getWorkingDir())
//    /* init */,workingDirLabel(new QLabel(QString::fromStdString(traitsIO.simulationFolder)))
    /* init */,bcActor(new BicrystalActor(renderWindow,renderer))
//    /* init */,mesh(traitsIO.meshFile,
//                    TextFileParser(traitsIO.polyFile).readMatrix<double>("A",3,3,true),
//                    TextFileParser(traitsIO.polyFile).readMatrix<double>("x0",1,3,true).transpose(),
//                    TextFileParser(traitsIO.polyFile).template readSet<int>("periodicFaceIDs",true))
//    /* init */,meshActor(new SimplicialMeshActor(renderWindow,renderer,mesh))
//    /* init */,ddConfigVtk(new DDconfigVtk(traitsIO.evlFolder,traitsIO.auxFolder,renderWindow,renderer,mesh))
    {
        renderer->SetBackground(1,1,1);
        renderWindow->AddRenderer(renderer);
        openglWidget->setRenderWindow(renderWindow);

        tabWidget->addTab(bcActor, tr(std::string("Bicrystal").c_str()));
//        tabWidget->addTab(meshActor, tr(std::string("Mesh").c_str()));

        mainLayout->addWidget(loadBicrystalButton,0,0,1,2);
        mainLayout->addWidget(tabWidget,1,0,1,1);
        mainLayout->addWidget(openglWidget,1,1,1,1);
        mainLayout->setColumnStretch(0, 3);
        mainLayout->setColumnStretch(1, 7);
        this->setLayout(mainLayout);
        
        connect(loadBicrystalButton,SIGNAL(released()), this, SLOT(getBicrystalFromFile()));

    }
    

    void oiViewerVTKwidget::getBicrystalFromFile()
    {
       const std::string fileName=QFileDialog::getOpenFileName(this, tr("Open File"),
                                                             "/Users/giacomo/Documents/oILAB/",
                                            tr("Images (*.txt)") ).toStdString();

        TextFileParser parser(fileName);
        const auto A(parser.readMatrix<double,3,3>("A",true));
//        if(A.rows()==3 && A.cols()==3)
//        {
            const auto R1(parser.readMatrix<double,3,3>("R1",true));
            const auto R2(parser.readMatrix<double,3,3>("R2",true));
            
            latA.reset(new Lattice<3>(A,R1));
            latB.reset(new Lattice<3>(A,R2));
            bc.reset(new BiCrystal<3>(*latA,*latB));
            bcActor->updateConfiguration(bc);
//        }
        
    
    }


} // namespace model

#endif







