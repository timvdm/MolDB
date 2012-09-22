#include "messagehandler.h"
#include "util.h"

#include <iostream>

namespace MolDB {

  std::string MessageHandler::type2string(enum MessageType type) const
  {
    switch (type) {
      case Information:
        return "Information";
      case Warning:
        return "Warning";
      case Error:
        return "Error";
    }

    return std::string();
  }

  MessageHandler::~MessageHandler()
  {
  }

  StdOutMessageHandler::StdOutMessageHandler(bool useColor) : m_color(useColor)
  {
  }

  void StdOutMessageHandler::report(enum MessageType type, const std::string &message)
  {
    if (m_color) {
      const char *red    = "\033[1;31m";
      const char *green  = "\033[1;32m";
      const char *yellow = "\033[1;33m";
      const char *blue   = "\033[1;34m";
      const char *normal = "\033[0m";

      switch (type) {
        case Warning:
          std::cout << yellow << type2string(type) << ": " << message << normal << std::endl;
          return;
        case Error:
          std::cout << red << type2string(type) << ": " << message << normal << std::endl;
          return;
        default:
          std::cout << type2string(type) << ": " << green << message << normal << std::endl;
          return;
      }
    }
    
    std::cout << type2string(type) << ": " << message << std::endl;
  }

  StdStreamMessageHandler::StdStreamMessageHandler(std::ostream *os) : m_os(os)
  {
  }

  void StdStreamMessageHandler::report(enum MessageType type, const std::string &message)
  {
    *m_os << type2string(type) << ": " << message << std::endl;
  }
 
  RecordingMessageHandler::RecordingMessageHandler(bool recordInformation, bool recordWarnings, 
      bool recordErrors) : m_messageHandler(0), m_counts(3, 0), 
      m_recordInformation(recordInformation), m_recordWarnings(recordWarnings),
      m_recordErrors(recordErrors)
  {
  }
      
  RecordingMessageHandler::RecordingMessageHandler(MessageHandler *messageHandler, 
      bool recordInformation, bool recordWarnings, bool recordErrors)
      : m_messageHandler(messageHandler), m_counts(3, 0), m_recordInformation(recordInformation),
      m_recordWarnings(recordWarnings), m_recordErrors(recordErrors)

  {
  }
   
  void RecordingMessageHandler::report(enum MessageType type, const std::string &message)
  {
    m_counts[type]++;
    switch (type) {
      case Information:
        if (m_recordInformation)
          m_information.push_back(Message(type, message));
        break;
      case Warning:
        if (m_recordWarnings)
          m_warnings.push_back(Message(type, message));
        break;
      case Error:
        if (m_recordErrors)
          m_errors.push_back(Message(type, message));
        break;
    }

    if (m_messageHandler)
      m_messageHandler->report(type, message);
  }

  std::size_t RecordingMessageHandler::numInformation() const
  {
    return m_counts[Information];
  }

  std::size_t RecordingMessageHandler::numWarnings() const
  {
    return m_counts[Warning];
  }

  std::size_t RecordingMessageHandler::numErrors() const
  {
    return m_counts[Error];
  }

  const std::vector<RecordingMessageHandler::Message>& RecordingMessageHandler::information() const
  {
    return m_information;
  }

  const std::vector<RecordingMessageHandler::Message>& RecordingMessageHandler::warnings() const
  {
    return m_warnings;
  }

  const std::vector<RecordingMessageHandler::Message>& RecordingMessageHandler::errors() const
  {
    return m_errors;
  }

}
