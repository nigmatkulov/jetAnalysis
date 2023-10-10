#ifndef GENTRACK_H
#define GENTRACK_H

// ROOT headers
#include "TObject.h"

//_________________
class GenTrack : public TObject {
  public:
    GenTrack();
    virtual ~GenTrack();

  private:

    ClassDef(GenTrack, 1)
};

#endif // #define GENTRACK_H