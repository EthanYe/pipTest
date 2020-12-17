#include "QPipe.h"
#include "hydroment.h"

#include <QPainter>
#include <QPen>
#include <QtMath>

//! [0]
QPipe::QPipe(Hydroment *startHydro, Hydroment *endHydro, QGraphicsItem *parent)
    : QGraphicsLineItem(parent), myStartHydro(startHydro), myEndHydro(endHydro)
{
    setFlag(QGraphicsItem::ItemIsSelectable, true);
    setPen(QPen(myColor, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
}
//! [0]

//! [1]
QRectF QPipe::boundingRect() const
{
    qreal extra = (pen().width() + 20) / 2.0;

    return QRectF(line().p1(), QSizeF(line().p2().x() - line().p1().x(),
                                      line().p2().y() - line().p1().y()))
        .normalized()
        .adjusted(-extra, -extra, extra, extra);
}
//! [1]

//! [2]
QPainterPath QPipe::shape() const
{
    QPainterPath path = QGraphicsLineItem::shape();
    path.addPolygon(QPipeHead);
    return path;
}
//! [2]

//! [3]
void QPipe::updatePosition()
{
    QLineF line(mapFromItem(myStartItem, 0, 0), mapFromItem(myEndItem, 0, 0));
    setLine(line);
}
//! [3]

//! [4]
void QPipe::paint(QPainter *painter, const QStyleOptionGraphicsItem *,
                  QWidget *)
{
    if (myStartItem->collidesWithItem(myEndItem))
        return;

    QPen myPen = pen();
    myPen.setColor(myColor);
    qreal QPipeSize = 20;
    painter->setPen(myPen);
    painter->setBrush(myColor);
//! [4] //! [5]

    QLineF centerLine(myStartItem->pos(), myEndItem->pos());
    QPolygonF endPolygon = myEndItem->polygon();
    QPointF p1 = endPolygon.first() + myEndItem->pos();
    QPointF intersectPoint;
    for (int i = 1; i < endPolygon.count(); ++i) {
        QPointF p2 = endPolygon.at(i) + myEndItem->pos();
        QLineF polyLine = QLineF(p1, p2);
        QLineF::IntersectionType intersectionType =
            polyLine.intersects(centerLine, &intersectPoint);
        if (intersectionType == QLineF::BoundedIntersection)
            break;
        p1 = p2;
    }

    setLine(QLineF(intersectPoint, myStartItem->pos()));
//! [5] //! [6]

    double angle = std::atan2(-line().dy(), line().dx());

    QPointF QPipeP1 = line().p1() + QPointF(sin(angle + M_PI / 3) * QPipeSize,
                                    cos(angle + M_PI / 3) * QPipeSize);
    QPointF QPipeP2 = line().p1() + QPointF(sin(angle + M_PI - M_PI / 3) * QPipeSize,
                                    cos(angle + M_PI - M_PI / 3) * QPipeSize);

    QPipeHead.clear();
    QPipeHead << line().p1() << QPipeP1 << QPipeP2;
//! [6] //! [7]
    painter->drawLine(line());
    painter->drawPolygon(QPipeHead);
    if (isSelected()) {
        painter->setPen(QPen(myColor, 1, Qt::DashLine));
        QLineF myLine = line();
        myLine.translate(0, 4.0);
        painter->drawLine(myLine);
        myLine.translate(0,-8.0);
        painter->drawLine(myLine);
    }
}
//! [7]
