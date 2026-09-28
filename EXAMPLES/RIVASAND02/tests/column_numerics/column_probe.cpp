// Optional diagnostic Tcl extension. Calls existing Domain/Element methods;
// it neither replaces the material implementation nor advances a solver.
#include <tcl.h>
#include <elementAPI.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <cstring>

static int probe(ClientData, Tcl_Interp* interp, int argc, const char* argv[]) {
    Domain* domain = OPS_GetDomain();
    if (!domain || argc < 2) return TCL_ERROR;
    int code = 0;
    if (argc == 3 && std::strcmp(argv[1], "materialIdentity") == 0) {
        int tag;
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) return TCL_ERROR;
        NDMaterial* material = OPS_getNDMaterial(tag);
        if (!material) return TCL_ERROR;
        NDMaterial* copy = material->getCopy("ThreeDimensional");
        if (!copy) return TCL_ERROR;
        const bool valid = std::strcmp(material->getClassType(), copy->getClassType()) == 0
            && material->getClassTag() == copy->getClassTag();
        Tcl_SetObjResult(interp, Tcl_NewStringObj(copy->getClassType(), -1));
        delete copy;
        return valid ? TCL_OK : TCL_ERROR;
    }
    if (std::strcmp(argv[1], "update") == 0) code = domain->update();
    else if (std::strcmp(argv[1], "commit") == 0) code = domain->commit();
    else if (std::strcmp(argv[1], "revert") == 0) code = domain->revertToLastCommit();
    else if (argc == 4 && std::strcmp(argv[1], "matrix") == 0) {
        int tag;
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) return TCL_ERROR;
        Element* element = domain->getElement(tag);
        if (!element) return TCL_ERROR;
        const Matrix* matrix = nullptr;
        if (std::strcmp(argv[3], "damping") == 0) matrix = &element->getDamp();
        else if (std::strcmp(argv[3], "tangent") == 0) matrix = &element->getTangentStiff();
        else if (std::strcmp(argv[3], "mass") == 0) matrix = &element->getMass();
        else return TCL_ERROR;
        Tcl_Obj* result = Tcl_NewListObj(0, nullptr);
        for (int i=0; i<matrix->noRows(); ++i)
            for (int j=0; j<matrix->noCols(); ++j)
                Tcl_ListObjAppendElement(interp, result, Tcl_NewDoubleObj((*matrix)(i,j)));
        Tcl_SetObjResult(interp, result);
        return TCL_OK;
    } else return TCL_ERROR;
    Tcl_SetObjResult(interp, Tcl_NewIntObj(code));
    return code == 0 ? TCL_OK : TCL_ERROR;
}

extern "C" int Columnprobe_Init(Tcl_Interp* interp) {
    Tcl_CreateCommand(interp, "columnProbe", probe, nullptr, nullptr);
    return Tcl_PkgProvide(interp, "columnProbe", "1.0");
}
