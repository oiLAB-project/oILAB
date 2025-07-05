//
// Created by Nikhil Chandra Admal on 7/3/25.
//

#ifndef OILAB_PYLATTICEMODULE_H
#define OILAB_PYLATTICEMODULE_H

namespace pyoilab
{
    template <int dim>
    class PyLatticeVector;

    template <int dim>
    class PyLatticeDirection;

    template <int dim>
    class PyReciprocalLatticeVector;

    template <int dim>
    class PyReciprocalLatticeDirection;
}

#include <LatticeBindings.h>
#include <LatticeVectorBindings.h>
#include <LatticeDirectionBindings.h>
#include <ReciprocalLatticeVectorBindings.h>
#include <ReciprocalLatticeDirectionBindings.h>

#endif //OILAB_PYLATTICEMODULE_H
