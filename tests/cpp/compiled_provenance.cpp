#include "tpsf/compiled_provenance.hpp"
#include <cassert>
#include <string>
int main() {
    const auto p = tpsf::compiled_engineering_provenance();
    assert(!p.statefinder_revision.empty());
    assert(!p.solver_revision.empty());
    assert(!p.compiler.empty());
    assert(!p.build_type.empty());
    return 0;
}
