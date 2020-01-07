/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2015 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
   */

#ifndef MODULES_FIELDS_EDITMESHBOUNDINGBOX_H
#define MODULES_FIELDS_EDITMESHBOUNDINGBOX_H

#include "Graphics/Widgets/BoundingBoxWidget.h"
#include <Dataflow/Network/GeometryGeneratingModule.h>
#include <Graphics/Datatypes/GeometryImpl.h>
#include <Graphics/Widgets/Widget.h>
#include <Core/Datatypes/Geometry.h>
#include <Modules/Fields/share.h>

namespace SCIRun {
  namespace Modules {
    namespace Fields {

      class EditMeshBoundingBoxImpl;

      class SCISHARE EditMeshBoundingBox : public Dataflow::Networks::GeometryGeneratingModule,
        public Has1InputPort<FieldPortTag>,
        public Has3OutputPorts < FieldPortTag, GeometryPortTag, MatrixPortTag >
      {
      public:
        static const int mDIMENSIONS = 3;
        EditMeshBoundingBox();
        void execute() override;
        void setStateDefaults() override;

        static const Core::Algorithms::AlgorithmParameterName TransformMatrix;
        static const Core::Algorithms::AlgorithmParameterName FieldTransformMatrix;
        static const Core::Algorithms::AlgorithmParameterName ResetToInput;
        static const Core::Algorithms::AlgorithmParameterName ResetSize;
        static const Core::Algorithms::AlgorithmParameterName ResetCenter;
        static const Core::Algorithms::AlgorithmParameterName DataSaved;
        //Input Field Attributes
        static const Core::Algorithms::AlgorithmParameterName InputCenterX;
        static const Core::Algorithms::AlgorithmParameterName InputCenterY;
        static const Core::Algorithms::AlgorithmParameterName InputCenterZ;
        static const Core::Algorithms::AlgorithmParameterName InputSizeX;
        static const Core::Algorithms::AlgorithmParameterName InputSizeY;
        static const Core::Algorithms::AlgorithmParameterName InputSizeZ;
        //Output Field Atributes
        static const Core::Algorithms::AlgorithmParameterName SetOutputCenter;
        static const Core::Algorithms::AlgorithmParameterName SetOutputSize;
        static const Core::Algorithms::AlgorithmParameterName OutputCenterX;
        static const Core::Algorithms::AlgorithmParameterName OutputCenterY;
        static const Core::Algorithms::AlgorithmParameterName OutputCenterZ;
        static const Core::Algorithms::AlgorithmParameterName OutputSizeX;
        static const Core::Algorithms::AlgorithmParameterName OutputSizeY;
        static const Core::Algorithms::AlgorithmParameterName OutputSizeZ;
        //Widget Scale/Mode
        static const Core::Algorithms::AlgorithmParameterName Scale;
        static const Core::Algorithms::AlgorithmParameterName OldScale;
        static const Core::Algorithms::AlgorithmParameterName ScaleChanged;
        static const Core::Algorithms::AlgorithmParameterName NoTranslation;
        static const Core::Algorithms::AlgorithmParameterName XYZTranslation;
        static const Core::Algorithms::AlgorithmParameterName RDITranslation;
        static const Core::Algorithms::AlgorithmParameterName RestrictX;
        static const Core::Algorithms::AlgorithmParameterName RestrictY;
        static const Core::Algorithms::AlgorithmParameterName RestrictZ;
        static const Core::Algorithms::AlgorithmParameterName RestrictR;
        static const Core::Algorithms::AlgorithmParameterName RestrictD;
        static const Core::Algorithms::AlgorithmParameterName RestrictI;

        static const Core::Algorithms::AlgorithmParameterName BoxMode;
        static const Core::Algorithms::AlgorithmParameterName BoxRealScale;

        INPUT_PORT(0, InputField, Field);
        OUTPUT_PORT(0, OutputField, Field);
        OUTPUT_PORT(1, Transformation_Widget, GeometryObject);
        OUTPUT_PORT(2, Transformation_Matrix, Matrix);

        MODULE_TRAITS_AND_INFO(ModuleHasUI)

      private:
        void clearVals();
        void computeWidgetBox(const Core::Geometry::BBox& box);
        void buildGeometryObject();
        void updateState(FieldHandle field);
        void sendOutputPorts();
        void resetToInputField();
        void processWidgetFeedback(const Core::Datatypes::ModuleFeedback& var);
        void adjustGeometryFromTransform(Core::Datatypes::MouseButton btn,
                                         Core::Datatypes::WidgetMovement move,
                                         const Core::Geometry::Transform& feedbackTrans);
        void generateGeomsList();
        void saveToParameters();
        void loadFromParameters();
        std::string convertForLabel(double coord);
        void updateInputFieldAttributes();

        Core::Geometry::Point pos_;
        std::vector<Core::Geometry::Vector> eigvecs_;
        std::vector<double> eigvals_;
        FieldHandle outputField_;

        std::vector<Graphics::Datatypes::GeometryHandle> geoms_;
        Core::Geometry::Transform trans_;
        Core::Geometry::Transform widgetAxes_;
        Core::Geometry::Transform widgetAxesOrthonormal_;
        Core::Geometry::Transform fieldTrans_;
        bool widgetMoved_;
        bool widgetAxesRotated_;
        bool firstRun_;
      };
    }
  }
}

#endif
