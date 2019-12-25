/*
For more information, please see: http://software.sci.utah.edu

The MIT License

Copyright (c) 2015 Scientific Computing and Imaging Institute,
University of Utah.


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

#ifndef Graphics_Widgets_BoundingBoxWidget_H
#define Graphics_Widgets_BoundingBoxWidget_H

#include <Core/Datatypes/Legacy/Field/FieldFwd.h>
#include <Core/GeometryPrimitives/GeomFwd.h>
#include <Graphics/Widgets/Widget.h>
#include <Graphics/Widgets/share.h>

namespace SCIRun {
  namespace Graphics {
    namespace Datatypes {
      enum BoundingBoxWidgetSection {
        BOX,
        CORNER_SCALE,
        FACE_ROTATE,
        FACE_SCALE,
        X_PLUS,
        Y_PLUS,
        Z_PLUS,
        X_MINUS,
        Y_MINUS,
        Z_MINUS };

      class SCISHARE BoundingBoxWidget : public CompositeWidget
      {
      public:
        BoundingBoxWidget(const Core::GeometryIDGenerator& idGenerator, const std::string& name,
                          double scale, const Core::Geometry::Point& center,
                          const std::vector<Core::Geometry::Vector>& eigvecs_,
                          const std::vector<double>& eigvals_,
                          int widget_num, int widget_iter);
        BoundingBoxWidget(const Core::GeometryIDGenerator& idGenerator, const std::string& name,
                          double scale, const BoxPosition& pos, const Core::Geometry::Point& center,
                          int widget_num, int widget_iter, const Core::Geometry::BBox& bbox);

        BoundingBoxWidget(const Core::GeometryIDGenerator& idGenerator, const std::string& name,
                          double scale, const Core::Geometry::Transform& trans,
                          const Core::Geometry::Point& center, int widget_num, int widget_iter);

      private:
        std::vector<WidgetHandle> translateWidgets_;
        std::vector<WidgetHandle> rotateWidgets_;
        std::vector<WidgetHandle> scaleWidgets_;
        std::vector<std::vector<std::vector<WidgetHandle>>> scaleAxisWidgets_;
        std::vector<std::string> allIds_;
        std::vector<std::string> translateIds_;
        std::vector<std::vector<std::vector<std::string>>> translateIdsByFace_;
        std::vector<std::vector<std::string> > translateIdsBySide_;
        std::vector<std::vector<std::vector<std::string>>> rotateIdsByFace_;
        std::vector<std::vector<std::vector<std::string>>> scaleIdsByFace_;
        std::vector<std::string> rotateIds_;
        std::vector<std::string> scaleIds_;
        std::vector<std::vector<std::vector<std::string>>> scaleAxisIds_;
        std::vector<std::vector<std::pair<WidgetMovement, std::vector<std::string>>>> translateMaps_;
        std::vector<std::pair<WidgetMovement, std::vector<std::string>>> rotateMap_;
        std::vector<std::pair<WidgetMovement, std::vector<std::string>>> scaleMap_;
        std::vector<std::vector<std::vector<std::pair<WidgetMovement, std::vector<std::string>>>>> scaleAxisMaps_;
        std::vector<std::vector<std::vector<std::pair<WidgetMovement, std::vector<std::string>>>>> scaleAxisUnidirectionalMaps_;
        const int DIMENSIONS_ = 3;
        const int EDGES_ = 12;
        const int CORNERS_ = 8;
        const int FACES_ = 6;
        std::string widgetName(size_t i, size_t id, size_t iter);

        int resolution_ = 10;
        const std::string deflPointCol_ {Core::Datatypes::ColorRGB(0.54, 0.6, 1.0).toString()};
        const std::string deflCol_ {Core::Datatypes::ColorRGB(0.5, 0.5, 0.5).toString()};
        const std::string resizeCol_ {Core::Datatypes::ColorRGB(0.54, 1.0, 0.60).toString()};
        std::string diskCol_ = resizeCol_;

        int widgetsIndex_;
        double scale_;
        double smallestEigval_;
        Core::Geometry::BBox bbox_;
        Core::Geometry::Point center_;
        std::vector<double> eigvals_;
        std::vector<Core::Geometry::Vector> eigvecs_;
        std::vector<Core::Geometry::Vector> scaledEigvecs_;
        std::vector<Core::Geometry::Point> corners_;
        std::vector<Core::Geometry::Point> facesStart_;
        std::vector<Core::Geometry::Point> facesEnd_;
        const float boxScale_ {1.0};
        const float rotSphereScale_ {2.0};
        const float resizeSphereScale_ {1.5};
        const float diskWidth_ {0.3};
        const float diskRadius_ {1.0};

        void addBox(const Core::GeometryIDGenerator& idGenerator, int widgetNum, int widgetIter);
        void createWidgets(const Core::GeometryIDGenerator& idGenerator, int widgetNum,
                           int widgetIter);
        void initWidgetCreation(const Core::GeometryIDGenerator& idGenerator, int widgetNum,
                                int widgetIter);
        void getEigenValuesAndEigenVectors();
        void getCorners();
        void getFacesStart();
        void getFacesEnd();
        void addCornerSpheres(const Core::GeometryIDGenerator& idGenerator,
                              int widgetNum, int widgetIter);
        void addFaceSphere(const Core::GeometryIDGenerator& idGenerator, int widgetNum,
                           int widgetIter);
        void addFaceCylinder(const Core::GeometryIDGenerator& idGenerator,
                             glm::mat4& scaleTrans, int widgetNum, int widgetIter);
        void getTranslateIds();
        void getRotateIds();
        void getScaleIds();
        void getScaleAxisIds();
      };

      using BoundingBoxWidgetHandle = SharedPointer<BoundingBoxWidget>;
    }
  }
}
#endif
