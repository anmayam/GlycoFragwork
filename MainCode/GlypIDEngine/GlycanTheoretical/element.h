// Written by Yin Wu (wuyin@indiana.edu)
// -------------------------------------------------------------------------------
// 
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at 
// http://www.apache.org/licenses/LICENSE-2.0

// modified by Anoop for GlypID


#ifndef ELEMENT_H
#define ELEMENT_H


/*
#include <qcolor.h>
#include <qnamespace.h>
#include <qstring.h>
#include <qvaluevector.h>
#include <system.h>

*/


namespace Engine
{
	namespace GlycanTheoretical
	{
		/*

		class Element;

		typedef QValueVector<Element> ElementVector;
		typedef list<Element> ElementList;

		//  Elements are valid if they have a value which is > EPSILON.

		class Element
		{

			private:
					void init( double value, QColor valueColor, int valuePattern, const QString& label, QColor labelColor );
					double m_value;
					QColor m_valueColor;
					int m_valuePattern;
					QString m_label;
					QColor m_labelColor;
					double m_propoints[2 * MAX_PROPOINTS];
		
			public:
				enum { INVALID = -1 };
				enum { NO_PROPORTION = -1 };
				enum { MAX_PROPOINTS = 3 }; // One proportional point per chart type

				Element( double value = INVALID, QColor valueColor = Qt::blue, int valuePattern = Qt::SolidPattern, const QString& label = QString::null, QColor labelColor = Qt::black) 
				{
							init( value, valueColor, valuePattern, label, labelColor );
							for ( int i = 0; i < MAX_PROPOINTS * 2; ++i )
								m_propoints[i] = NO_PROPORTION;
				}
				
				~Element() {}

				bool isValid() const { return m_value > EPSILON; }

				double value() const { return m_value; }
				QColor valueColor() const { return m_valueColor; }
				int valuePattern() const { return m_valuePattern; }
				QString label() const { return m_label; }
				QColor labelColor() const { return m_labelColor; }
				double proX( int index ) const;
				double proY( int index ) const;

				void set( double value = INVALID, QColor valueColor = Qt::gray, int valuePattern = Qt::SolidPattern, const QString& label = QString::null, QColor labelColor = Qt::black)
				{
					init( value, valueColor, valuePattern, label, labelColor );
				}
				void setValue( double value ) { m_value = value; }
				void setValueColor( QColor valueColor ) { m_valueColor = valueColor; }
				void setValuePattern( int valuePattern );
				void setLabel( const QString& label ) { m_label = label; }
				void setLabelColor( QColor labelColor ) { m_labelColor = labelColor; }
				void setProX( int index, double value );
				void setProY( int index, double value );

				#ifdef Q_FULL_TEMPLATE_INSTANTIATION
					// xlC 3.x workaround
					Q_DUMMY_COMPARISON_OPERATOR(Element)
					bool operator!=( const Element& e) const {
						return ( !(e == *this) );
					}
				#endif

				double x_value, y_value;
				string message;

				
		};


		QTextStream &operator<<( QTextStream&, const Element& );
		QTextStream &operator>>( QTextStream&, Element& );

		*/
	}
}

#endif
