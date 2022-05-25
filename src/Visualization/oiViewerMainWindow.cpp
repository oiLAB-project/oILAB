/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_oiViewerMainWindow_cpp_
#define model_oiViewerMainWindow_cpp_

#include <iostream>

#include <oiViewerMainWindow.h>
#include <QMenu>
#include <QMenuBar>


namespace gbLAB
{
    
    oiViewerMainWindow::oiViewerMainWindow() :
    /* init */ viewerCount(0)
    /* init */,tabWidget(new QTabWidget(this))
    /* init */,newViewerAction(new QAction(tr("&New"), this))
    {
        resize(1920, 1080);

        QMenu* viewertMenu=menuBar()->addMenu(tr("&Viewer"));
        viewertMenu->addAction(newViewerAction);
        connect(newViewerAction, &QAction::triggered, this, &oiViewerMainWindow::newViewerTab);
        
        setCentralWidget(tabWidget);
    }

    void oiViewerMainWindow::newViewerTab()
    {
        try
        {
//            tabWidget->addTab(new oiViewerVTKwidget(this), tr(std::string("Viewer "+std::to_string(viewerCount)).c_str()));
            viewerCount++;
        }
        catch(const std::exception& e)
        {
            std::cout<<e.what()<<std::endl;
        }
    }

} // namespace model
#endif







