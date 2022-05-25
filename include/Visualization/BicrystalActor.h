/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BicrystalActor_h_
#define model_BicrystalActor_h_


#include <deque>
#include <string>
#include <memory>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QLineEdit>

#include <vtkGenericOpenGLRenderWindow.h>
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

#include <BiCrystal.h>


namespace gbLAB
{
    struct BicrystalActor : public QWidget
//: public DDconfigVtkBase
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        vtkRenderer* const renderer;

        QGridLayout* mainLayout;
        QCheckBox* showA;
        QCheckBox* showB;

//        std::shared_ptr<Lattice<3>> latticeA;
//        std::shared_ptr<Lattice<3>> latticeB;
        
//        QCheckBox* showNodeLabels;
//        QCheckBox* showVelocities;
//        QLineEdit* velocityScaleEdit;

        
        public:
                
        vtkSmartPointer<vtkPolyData> aPolyData;
        vtkSmartPointer<vtkGlyph3D> aGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> aMapper;
        vtkSmartPointer<vtkActor> aActor;

        vtkSmartPointer<vtkPolyData> bPolyData;
        vtkSmartPointer<vtkGlyph3D> bGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> bMapper;
        vtkSmartPointer<vtkActor> bActor;

        
//        vtkSmartPointer<vtkPolyData> labelPolyData;
//        vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
//        vtkSmartPointer<vtkActor2D> labelActor;
//
//        vtkSmartPointer<vtkPolyData> velocityPolyData;
//        vtkSmartPointer<vtkGlyph3D> velocityGlyphs;
//        vtkSmartPointer<vtkPolyDataMapper> velocityMapper;
//        vtkSmartPointer<vtkActor> velocityActor;
//
//        vtkSmartPointer<vtkPolyData> singleNodeLabelPolyData;
//        vtkSmartPointer<vtkLabeledDataMapper> singleNodeLabelMapper;
//        vtkSmartPointer<vtkActor2D> singleNodeLabelActor;
//
//        size_t singleNodeID;
//        unsigned char nodeClr[4][3];
        
        BicrystalActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const);
        void updateConfiguration(const std::shared_ptr<BiCrystal<3>>& bc);
        
        
    };
    
} // namespace model
#endif
