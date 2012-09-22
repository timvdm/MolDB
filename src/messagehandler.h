#ifndef MOLDB_MESSAGEHANDLER_H
#define MOLDB_MESSAGEHANDLER_H

#include <string>
#include <iosfwd>
#include <vector>

namespace MolDB {

  enum MessageType {
    Information,
    Warning,
    Error
  };

  class MessageHandler
  {
    public:
      std::string type2string(enum MessageType type) const;

      virtual void report(enum MessageType type, const std::string &message) = 0;
      virtual ~MessageHandler();
  };

  class StdOutMessageHandler : public MessageHandler
  {
    public:
      StdOutMessageHandler(bool useColor = true);
      void report(enum MessageType type, const std::string &message);
    private:
      bool m_color;
  };

  class StdStreamMessageHandler : public MessageHandler
  {
    public:
      StdStreamMessageHandler(std::ostream *os);
      void report(enum MessageType type, const std::string &message);
    private:
      std::ostream *m_os;
  };

  class RecordingMessageHandler : public MessageHandler
  {
    public:
      struct Message
      {
        Message(enum MessageType type_, const std::string &message_)
            : type(type_), message(message_)
        {
        }

        enum MessageType type;
        std::string message;
      };

      RecordingMessageHandler(bool recordInformation = false, bool recordWarnings = true, bool recordErrors = true);
      RecordingMessageHandler(MessageHandler *messageHandler, bool recordInformation = false, bool recordWarnings = true, bool recordErrors = true);
      void report(enum MessageType type, const std::string &message);

      std::size_t numInformation() const;
      std::size_t numWarnings() const;
      std::size_t numErrors() const;

      const std::vector<Message>& information() const;
      const std::vector<Message>& warnings() const;
      const std::vector<Message>& errors() const;

    private:
      MessageHandler *m_messageHandler;
      std::vector<std::size_t> m_counts;
      std::vector<Message> m_information;
      std::vector<Message> m_warnings;
      std::vector<Message> m_errors;
      bool m_recordInformation;
      bool m_recordWarnings;
      bool m_recordErrors;
  };


}

#endif
