/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_oiViewerVTKwidget_h_
#define model_oiViewerVTKwidget_h_

#include <QGridLayout>
#include <QLabel>
#include <QTabWidget>
#include <QPushButton>

#include <QVTKOpenGLStereoWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>

#include <BiCrystal.h>
#include <TextFileParser.h>
#include <BicrystalActor.h>

namespace gbLAB
{
    

    

struct oiViewerVTKwidget : public QWidget
//public QVTKOpenGLStereoWidget
{
    Q_OBJECT
    

    
private:
    
    QGridLayout* mainLayout;
    QPushButton* loadBicrystalButton;
    QTabWidget* tabWidget;
    QVTKOpenGLStereoWidget* openglWidget;
    
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderer> renderer;

    std::shared_ptr<Lattice<3>> latA;
    std::shared_ptr<Lattice<3>> latB;
    std::shared_ptr<BiCrystal<3>> bc;
    BicrystalActor* bcActor;

    
    private slots:
    void getBicrystalFromFile();
    

    
public:
    
//    const std::string workingDir;
//    const DDtraitsIO traitsIO;
//    QLabel* workingDirLabel;
//    Lattice<3> L1(A);
//    Lattice<3> L2(A,R);
//    BiCrystal<3> bc;
    
    
    oiViewerVTKwidget(QWidget *parent);
    

};


} // namespace model

#endif







