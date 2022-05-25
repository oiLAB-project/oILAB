/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_oiViewerMainWindow_H_
#define model_oiViewerMainWindow_H_

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QWidget>
#include <QGridLayout>
#include <QPushButton>
#include <QTabWidget>
#include <QAction>


#include <oiViewerVTKwidget.h>

namespace gbLAB
{
    
struct oiViewerMainWindow : public QMainWindow
{
    // https://kitware.github.io/vtk-examples/site/Cxx/Qt/RenderWindowNoUiFile/
    // https://kitware.github.io/vtk-examples/site/Cxx/Qt/SideBySideRenderWindowsQt/
    //https://vtk.org/doc/nightly/html/classQVTKOpenGLNativeWidget.html

    Q_OBJECT
    
    private slots:
    
        void newViewerTab();
    
    
    private:
    
        size_t viewerCount;
        QTabWidget* tabWidget;
        QAction* newViewerAction;

    public:

        oiViewerMainWindow();

};

} // namespace model
#endif







